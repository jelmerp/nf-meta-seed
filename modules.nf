// QC reads
process FASTQC {
    publishDir "${params.outdir}/read-qc/fastqc", mode: 'copy', pattern: '*html'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc --quiet ${reads[0]} ${reads[1]}
    """
}

// QC and trim reads
process FASTP {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*fastq.gz'), emit: fastq
    path '*fastp_report*', emit: report

    script:
    """
    fastp \
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${sample_id}_R1_trimmed.fastq.gz \\
        -O ${sample_id}_R2_trimmed.fastq.gz \\
        -h ${sample_id}_fastp_report.html \\
        -j ${sample_id}_fastp_report.json \\
        -w ${task.cpus}
    """
}

// Count reads in FASTQ files
process READCOUNT {
    input:
    path reads
    val run_id

    output:
    path 'readcounts*tsv'

    script:
    """
    seqkit stats \
        --threads $task.cpus \
        --tabular \
        ${reads} \
        > readcounts_${run_id}.tsv
    """
}

process REPORT_READLOSS {
    publishDir "${params.outdir}/hostremove/readcounts", mode: 'copy', pattern: '*tsv'

    input:
    path readcount_pre
    path readcount_post

    output:
    path "readloss.tsv"

    script:
    """
    echo -e "sample\tpre_host_removal\tpost_host_removal\tpercent_removed" > readloss.tsv
    join \
        <(tail -n+2 ${readcount_pre} | cut -f 1,4 | sort -k1,1 | sed 's/_.*gz//') \
        <(tail -n+2 ${readcount_post} | cut -f 1,4 | sort -k1,1 | sed 's/_.*gz//') |
        awk '{print \$0, "\t", ((\$2-\$3)/\$2)*100}' \
        >> readloss.tsv
    """
}

// Download and extract a pre-made Kraken DB from a URL
process KRAKENDB_DL_DB {
    input:
    val db_url

    output:
    path 'kraken_db'

    when:
    params.skip_kraken == false

    script:
    """
    wget ${db_url}
    tar -xzvf *.tar.gz -C kraken_db
    """
}

// Download the NCBI taxonomy as a starting point for a custom Kraken DB
process KRAKENDB_DL_TAX {
    output:
    path 'kraken_db_tax'

    when:
    params.skip_kraken == false

    script:
    """
    kraken2-build --download-taxonomy --db kraken_db_tax
    """
}

// Download taxon-specific libraries as a starting point for a custom Kraken DB
process KRAKENDB_DL_LIB {
    input:
    val library

    output:
    path '*'

    when:
    params.skip_kraken == false

    script:
    """
    kraken2-build --download-library ${library} --db db_dir

    mv db_dir/library/* .
    rm -r db_dir
    """
}

// Combine separate dirs with libraries into one
process KRAKENDB_COMBINE_LIBS {
    input:
    path lib_dirs

    output:
    path 'lib_dir'

    when:
    params.skip_kraken == false

    script:
    """
    mkdir -p lib_dir
    cp -rLv ${lib_dirs} lib_dir/
    """
}

// Add a set of custom genomes in a dir to a custom Kraken DB
process KRAKENDB_COMBINE_AND_ADD {
    input:
    path tax_dir
    path lib_dir
    path genomes_to_add_dir

    output:
    path "kraken_db_unbuilt"

    when:
    params.skip_kraken == false

    script:
    """
    mkdir -p kraken_db_unbuilt/library
    cp -rLv ${tax_dir}/taxonomy kraken_db_unbuilt
    cp -rLv ${lib_dir}/* kraken_db_unbuilt/library/

    if [[ -d ${genomes_to_add_dir} ]]; then
        for fasta in ${genomes_to_add_dir}/*fna; do
            kraken2-build --add-to-library \$fasta --db kraken_db_unbuilt
        done
    fi
    """
}

// Build the final custom Kraken DB
process KRAKENDB_BUILD {
    input:
    path kraken_db_unbuilt

    output:
    path 'kraken_db'

    when:
    params.skip_kraken == false

    script:
    """
    mkdir -p kraken_db
    cp -rLv ${kraken_db_unbuilt}/* kraken_db/
    
    kraken2-build --build --db kraken_db --threads ${task.cpus}
    """
}

// Run Kraken to assign taxonomy to reads
process KRAKEN {
    input:
    tuple val(sample_id), path(reads)
    path kraken_db
    val confidence
    val minhitgroups
    val run_id

    output:
    tuple val(sample_id), path('*report.txt'), path('*main.txt'), emit: k_extract
    tuple val(sample_id), path('*report.txt'), emit: report
    tuple val(sample_id), path('*main.txt'), emit: main_out
    tuple val(sample_id), path('*fastq'), emit: classified_fq
    path 'mqc/*.txt', emit: mqc // For MultiQC, to avoid complications with tuple
    path '*.log'

    when:
    params.skip_kraken == false

    script:
    """
    kraken2 \\
        --db ${kraken_db} \\
        --report ${sample_id}_${run_id}_kraken-report.txt \\
        --output ${sample_id}_${run_id}_kraken-main.txt \\
        --classified-out ${sample_id}#.fastq \\
        --minimum-hit-groups ${minhitgroups} \\
        --confidence ${confidence} \\
        --gzip-compressed \\
        --paired \\
        --threads ${task.cpus} \\
        ${reads[0]} \\
        ${reads[1]}
    
    mkdir -p mqc
    cp *kraken-report.txt mqc/

    cp .command.log command_kraken_${sample_id}.log
    """
}

// Extract reads from Kraken run
process KRAKEN_EXTRACT {
    input:
    tuple val(sample_id), path(kraken_report), path(kraken_output), path(reads)
    val tax_ids

    output:
    tuple val(sample_id), path('*fastq.gz'), emit: fastq
    path 'logs'

    script:
    """
    extract_kraken_reads.py \\
        -t ${tax_ids} \\
        -k ${kraken_output} \\
        -r ${kraken_report} \\
        -s ${reads[0]} \\
        -s2 ${reads[1]} \\
        -o ${sample_id}_R1.fastq \\
        -o2 ${sample_id}_R2.fastq \\
        --exclude \\
        --include-children \\
        --fastq-output

    gzip -fv ${sample_id}_R1.fastq ${sample_id}_R2.fastq

    mkdir -p logs
    grep "reads printed to file" .command.log > nreads_extracted_${sample_id}.log
    grep -v "reads processed" .command.log > command_extract_kraken_${sample_id}.log
    """
}

// Build a Bracken DB
process BRACKENDB_BUILD {
    input:
    path kraken_db
    val read_len

    output:
    path 'bracken_db'

    when:
    params.skip_bracken == false

    script:
    """
    cp -rLv ${kraken_db} bracken_db

    bracken-build -d bracken_db -l ${read_len} -t ${task.cpus}
    """
}

// Get Krona Taxonomy
process KRONA_TAX {
    output:
    path 'tax.tab'

    script:
    """
    wget https://raw.githubusercontent.com/marbl/Krona/refs/heads/master/KronaTools/updateTaxonomy.sh
    wget https://raw.githubusercontent.com/marbl/Krona/refs/heads/master/KronaTools/scripts/taxonomy.make
    wget https://raw.githubusercontent.com/marbl/Krona/refs/heads/master/KronaTools/scripts/extractTaxonomy.pl
    mkdir -p scripts
    mv taxonomy.make extractTaxonomy.pl scripts/
    chmod +x scripts/*

    bash updateTaxonomy.sh tax.tab
    """
}

// Run Krona
process KRONA {
    input:
    //tuple val(sample_id), path(kraken_report), path(kraken_main)
    tuple val(sample_id), path(kraken_main)
    path taxfile

    output:
    path '*html'

    script:
    """
    ktImportTaxonomy \\
        -q 2 \\
        -t 3 \\
        -tax ${taxfile} \\
        ${kraken_main} \\
        -o krona_${sample_id}.html
    """
}

// Run Bracken
process BRACKEN {
    input:
    tuple val(sample_id), path(kraken_report)
    path bracken_db
    val tax_level
    val min_reads
    val read_len
    val run_id

    output:
    path '*bracken*out.txt', emit: main_out
    path '*bracken*report.txt', emit: report

    when:
    params.skip_bracken == false

    script:
    """
    bracken \
        -i ${kraken_report} \\
        -d ${bracken_db} \\
        -o ${sample_id}_${run_id}_bracken-out.txt \\
        -w ${sample_id}_${run_id}_bracken-report.txt \\
        -r ${read_len} \\
        -l ${tax_level} \\
        -t ${min_reads}
    """
}

process BIOM {
    input:
    path kraken_reports
    val run_id

    output:
    path "*biom"

    script:
    """
    kraken-biom --fmt json -o ${run_id}.biom ${kraken_reports}
    """
}

// Assemble reads
process ASSEMBLY {
    publishDir "${params.outdir}/assembly/spades/asm", mode: 'copy', pattern: '*scaffolds.fasta'
    publishDir "${params.outdir}/assembly/spades/logs", mode: 'copy', pattern: '*spades.log'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_scaffolds.fasta"), emit: fasta
    path "${sample_id}_contigs.fasta"
    path "${sample_id}_spades.log"

    script:
    def memory_gb = MemoryUnit.of("${task.memory}").toUnit('GB')
    """
    spades.py \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o outdir \\
        --only-assembler \\
        --meta \\
        --threads ${task.cpus} \\
        --memory ${memory_gb}
    
    mv outdir/contigs.fasta ${sample_id}_contigs.fasta
    mv outdir/scaffolds.fasta ${sample_id}_scaffolds.fasta
    mv outdir/spades.log ${sample_id}_spades.log
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path multiqc_input
    val run_id

    output:
    path 'multiqc*html'

    script:
    """
    multiqc \\
        --filename multiqc_${run_id}.html \\
        --interactive \\
        .
    """
}

process HOST_INDEX {
    publishDir "${params.outdir}/hostremove/host_index", mode: 'copy', enabled: params.save_host_index

    input:
    path host_fasta

    output:
    path 'host_index_dir'

    script:
    """
    bowtie2-build ${host_fasta} host_index

    mkdir -p host_index_dir
    mv *bt2 host_index_dir/
    """
}

process HOST_REMOVE_ALIGN {
    input:
    path host_index_dir
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*.fastq.gz'), emit: fastq
    path ('bowtie-log*txt'), emit: logs

    // Note: redirecting stdout to /dev/null : this is the mapped SAM, which is not needed

    script:
    """
    index_prefix=\$(ls ${host_index_dir} | head -n1 | sed -E "s/.[0-9]+.bt2//")
    index_prefix_full=${host_index_dir}/\$index_prefix

    bowtie2 \
        --threads ${task.cpus} \
        -x \$index_prefix_full \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        --local \
        --un-conc-gz ${sample_id}_host_removed_reads \
        > /dev/null \
        2> bowtie-log_${sample_id}.txt

    mv -v ${sample_id}_host_removed_reads.1 ${sample_id}_hostrm_R1.fastq.gz
    mv -v ${sample_id}_host_removed_reads.2 ${sample_id}_hostrm_R2.fastq.gz
    """
}

process MAXBIN2 {
    publishDir "${params.outdir}/assembly/maxbin2/bins", mode: "copy", pattern: "*.fasta"
    publishDir "${params.outdir}/assembly/maxbin2/cov_bin", mode: "copy", pattern: "*.abundance"

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path('*.fasta'), emit: fasta, optional: true
    path '*abundance', optional: true

    script:
    """
    set +e
    
    run_MaxBin.pl \\
        -contig ${assembly} \\
        -reads ${reads[0]} \\
        -reads2 ${reads[1]} \\
        -out ${sample_id}

    if [[ \$? -ne 0 ]]; then
        touch maxbin_dummy_${sample_id}.fasta
    fi
    """
}

process METABAT2 {
    publishDir "${params.outdir}/assembly/metabat2/cov_contig", mode: "copy", pattern: "*depth.txt"
    publishDir "${params.outdir}/assembly/metabat2/bins", mode: "copy", pattern: "*.fa"
    publishDir "${params.outdir}/assembly/metabat2/cov_bin", mode: "copy", pattern: "*BinInfo.txt"

    input:
    tuple val(sample_id), path(assembly), path(bam), path(bam_idx), path(bed)

    output:
    tuple val(sample_id), path('*.fa'), emit: fasta, optional: true
    path '*depth.txt', optional: true
    path '*BinInfo.txt', optional: true

    script:
    """
    jgi_summarize_bam_contig_depths \\
        --outputDepth ${sample_id}_depth.txt \\
        "${bam}"

    metabat2 \\
        -i "${assembly}" \\
        -a ${sample_id}_depth.txt \\
        -m 1500 \\
        --maxP 75 \\
        -s 100000 \\
        -o "${sample_id}"

    n_files=`find . -type f -name "*fa" | wc -l`
    if [[ \$n_files -eq 0 ]]; then
        touch metabat_dummy_${sample_id}.fa
    fi
    """
}

process CONCOCT {
    publishDir "${params.outdir}/assembly/concoct/cov", mode: "copy", pattern: "covtable_*"
    publishDir "${params.outdir}/assembly/concoct/bins", mode: "copy", pattern: "bins_*/*.fa"

    input:
    tuple val(sample_id), path(assembly), path(bam), path(bam_idx), path(bed)

    output:
    tuple val(sample_id), path('bins_*/*.fa'), emit: fasta, optional: true
    path 'covtable_*'
    // path 'contigs10k.fasta' 'clustering_gt1000.csv' 'merged.csv'

    //TODO: don't hardcode read length

    script:
    """
    cp "${bed}" copy.bed # This is needed so Nextflow won't modify the input via the symlink

    cut_up_fasta.py \\
        "${assembly}" \\
        -c 100000 \\
        --merge_last \\
        -b copy.bed \\
        > contigs10k_${sample_id}.fasta
    
    concoct_coverage_table.py \\
        copy.bed \\
        "${bam}" \\
        > covtable_${sample_id}.tsv
    
    concoct \\
        --composition_file contigs10k_${sample_id}.fasta \\
        --coverage_file covtable_${sample_id}.tsv \\
        --seed 100 \\
        --read_length 150 \\
        --length_threshold 1000 \\
        --threads ${task.cpus}
    
    merge_cutup_clustering.py \\
        clustering_gt1000.csv \\
        > merged_${sample_id}.csv

    mkdir -p bins_${sample_id}
    extract_fasta_bins.py \\
        "${assembly}" \\
        merged_${sample_id}.csv \\
        --output_path bins_${sample_id}
    """
}

process DREP {
    publishDir "${params.outdir}/assembly/drep", mode: "copy", pattern: "dereplicated_genomes"
    publishDir "${params.outdir}/assembly/drep/logs", mode: "copy", pattern: "*.log"

    input:
    //tuple val(sample_id), path(concoct_fa), path(maxbin_fa), path(metabat_fa)
    tuple val(sample_id), path(concoct_fa), path(metabat_fa)

    output:
    path 'drep_out/data_tables'
    path 'drep_out/dereplicated_genomes', emit: derepped
    path 'drep_out/data/checkM/checkM_outdir/results.tsv', emit: checkM_result

    when:
    params.skip_drep == false

    script:
    """
    dRep dereplicate drep_out -g *.f*a
    cp .command.log drep_${sample_id}.log
    """
}

process MAP2ASSEMBLY {
    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}.bed")

    script:
    """
    bwa index -p "${sample_id}" "${assembly}"

    bwa mem -t ${task.cpus} -a "${sample_id}" ${reads[0]} ${reads[1]} |
        samtools sort -o "${sample_id}".bam -

    samtools index "${sample_id}".bam

    bedtools bamtobed -i "${sample_id}".bam > "${sample_id}".bed
    """
}

process METAPHLAN_DB {
    output:
    path 'metaphlan_db'

    script:
    """
    metaphlan --install --bowtie2db metaphlan_db
    """
}

process METAPHLAN {
    publishDir "${params.outdir}/classif-read/metaphlan/by-sample", mode: "copy", pattern: "*_profile.txt"

    input:
    tuple val(sample_id), path(reads)
    path metaphlan_db

    output:
    tuple val(sample_id), path("*_profile.txt"), emit: profile
    tuple val(sample_id), path('*.bowtie2out.txt'), optional: true, emit: bt2out
    path "*_profile.txt", emit: mqc

    script:
    """
    DB_VERSION=`find -L ${metaphlan_db} -name "*.rev.1.bt2*" | xargs -I{} basename {} .rev.1.bt2l`
    echo \$DB_VERSION > metaphlan_db_version.txt

    metaphlan \\
        --input_type fastq \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --subsampling_paired 100000000 \\
        --db_dir ${metaphlan_db} \\
        --index \$DB_VERSION \\
        --mapout ${sample_id}_bowtie2out.txt \\
        --output_file ${sample_id}_profile.txt \\
        --nproc ${task.cpus}
    """
}

process METAPHLAN_MERGE {
    publishDir "${params.outdir}/classif-read/metaphlan/merged", mode: "copy", pattern: "*metaphlan.txt"

    input:
    path metaphlan_profiles

    output:
    path '*txt'

    script:
    """
    merge_metaphlan_tables.py \\
        -o merged_metaphlan.txt \\
        ${metaphlan_profiles}
    """
}

process SOURMASH_DB {
    output:
    path 'sourmash_db.zip'
    path 'sourmash_taxdb.csv'

    script:
    """
    SOURMASH_DB_URL=https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs220/gtdb-rs220-k31.dna.zip
    SOURMAX_TAXDB_URL=https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs220/gtdb-rs220.lineages.csv

    curl --insecure -JLsS -o sourmash_db.zip \$SOURMASH_DB_URL
    curl --insecure -JLsS -o sourmash_taxdb.csv \$SOURMAX_TAXDB_URL
    """
}

process SOURMASH {
    publishDir "${params.outdir}/classif-asm/sourmash", mode: "copy"

    input:
    tuple val(sample_id), path(assembly)
    path sourmash_db
    path sourmash_taxdb

    output:
    path '*txt'
    path '*csv'

    script:
    """
    KMER_SIZE=31
    THRESHOLD_BP=10000

    sourmash sketch dna \
        -p abund,k="\$KMER_SIZE" \
        --output ${sample_id}.sig \
        $assembly
    
    sourmash gather \
        --output ${sample_id}_gather.csv \
        --threshold-bp \$THRESHOLD_BP \
        ${sample_id}.sig
        $sourmash_taxdb

    sourmash tax metagenome \
        --gather-csv ${sample_id}_gather.csv \
        --taxonomy-csv $sourmash_taxdb \
        --use-abundances \
        --rank species \
        --output-format human \
        --output-format csv_summary \
        --output-format krona \
        --output-format lineage_summary \
        --output-format kreport \
        --output-dir . \
        --output-base $sample_id
    """
}
