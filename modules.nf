// QC reads
process FASTQC {
    publishDir "$params.outdir/fastqc", mode: 'copy', pattern: '*html'
    
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
        -w $task.cpus
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
    wget $db_url
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
    kraken2-build --download-library $library --db db_dir

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
    cp -rLv $lib_dirs lib_dir/
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
    cp -rLv $kraken_db_unbuilt/* kraken_db/
    
    kraken2-build --build --db kraken_db --threads $task.cpus
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
    tuple val(sample_id), path('*report.txt'), path('*main.txt'), emit: output
    tuple val(sample_id), path('*fastq'), emit: classified_fq
    path 'mqc/*.txt', emit: mqc // For MultiQC, to avoid complications with tuple
    path '*.log'

    when:
    params.skip_kraken == false

    script:
    """
    kraken2 \\
        --db ${kraken_db} \\
        --report ${sample_id}_kraken-report.txt \\
        --output ${sample_id}_kraken-main.txt \\
        --classified-out ${sample_id}#.fastq \\
        --minimum-hit-groups $minhitgroups \\
        --confidence $confidence \\
        --gzip-compressed \\
        --paired \\
        --threads $task.cpus \\
        ${reads[0]} \\
        ${reads[1]}
    
    mkdir -p mqc
    cp ${sample_id}_kraken-report.txt mqc/${run_id}_${sample_id}_kraken-report.txt

    cp .command.log command_kraken_${sample_id}.log
    """
}

// Extract reads from Kraken run
process KRAKEN_EXTRACT {
    input:
    tuple val(sample_id), path(kraken_report), path(kraken_output), path(reads)
    val tax_ids

    output:
    tuple val(sample_id), path('*fastq.gz'), emit: fq
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
    cp -rLv $kraken_db bracken_db

    bracken-build -d bracken_db -l $read_len -t $task.cpus
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
    tuple val(sample_id), path(kraken_report), path(kraken_main)
    path taxfile

    output:
    path '*html'

    //? [-q <integer>]   Column of input files to use as query ID. Required if magnitude files are specified. [Default: '1']
    //? [-t <integer>]   Column of input files to use as taxonomy ID. [Default: '2']
    
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
    publishDir "$params.outdir/bracken", mode: 'copy', pattern: '*txt'

    input:
    tuple val(sample_id), path(kraken_report), path(kraken_main)
    path bracken_db
    val tax_level
    val min_reads
    val read_len

    output:
    path '*bracken*txt'

    when:
    params.skip_bracken == false

    script:
    """
    bracken \
        -i ${kraken_report} \\
        -d ${bracken_db} \\
        -o ${sample_id}_bracken-out.txt \\
        -w ${sample_id}_bracken-report.txt \\
        -r ${read_len} \\
        -l ${tax_level} \\
        -t ${min_reads}
    """
}

// Assemble reads
process ASSEMBLY {
    publishDir "$params.outdir/spades", mode: 'copy', pattern: '*scaffolds.fasta'
    publishDir "$params.outdir/spades", mode: 'copy', pattern: '*spades.log'

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
        --threads $task.cpus \\
        --memory $memory_gb
    
    mv outdir/contigs.fasta ${sample_id}_contigs.fasta
    mv outdir/scaffolds.fasta ${sample_id}_scaffolds.fasta
    mv outdir/spades.log ${sample_id}_spades.log
    """
}

process MULTIQC {
    publishDir "$params.outdir/multiqc", mode: 'copy'

    input:
    path multiqc_input

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc --interactive .
    """
}

process HOST_INDEX {
    publishDir "${params.outdir}/hot_index", mode: 'copy', enabled: params.save_host_index

    input:
    path host_fasta

    output:
    path 'host_index_dir'

    script:
    """
    bowtie2-build $host_fasta host_index

    mkdir -p host_index_dir
    mv *bt2 host_index_dir/
    """
}

process HOST_REMOVE_ALIGN {
    input:
    path(host_index_dir)
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*.fastq.gz'), emit: fastq
    path('bowtie-log*txt')

    shell:
    '''
    index_prefix=$(ls !{host_index_dir} | head -n1 | sed -E "s/.[0-9]+.bt2//")
    index_prefix_full=!{host_index_dir}/$index_prefix

    bowtie2 \\
        -p !{task.cpus} \\
        -x $index_prefix_full \\
        -1 !{reads[0]} \\
        -2 !{reads[1]} \\
        --local \\
        --un-conc !{sample_id}_host_removed_reads \\
        -S !{sample_id}_mapped_unmapped.sam \\
        2> bowtie-log_!{sample_id}.txt

    gzip -c !{sample_id}_host_removed_reads.1 > !{sample_id}_hostrm_R1.fastq.gz
    gzip -c !{sample_id}_host_removed_reads.2 > !{sample_id}_hostrm_R2.fastq.gz
    '''
}

process MAXBIN2 {
    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path('*.fasta'), emit: fasta

    script:
    """
    set +e
    
    run_MaxBin.pl \\
        -contig $assembly \\
        -reads ${reads[0]} \\
        -reads2 ${reads[1]} \\
        -out $sample_id

    if [[ \$? -ne 0 ]]; then
        touch maxbin_dummy_${sample_id}.fasta
    fi
    """
}

process METABAT2 {
    input:
    tuple val(sample_id), path(assembly), path(bam), path(bam_idx), path(bed)

    output:
    tuple val(sample_id), path('*.fa'), emit: fasta, optional: true
    path 'depth.txt'

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt "$bam"

    metabat2 -i "$assembly" -a depth.txt -m 1500 --maxP 75 -s 100000 -o "$sample_id"

    n_files=`find . -type f -name "*fa" | wc -l`
    if [[ \$n_files -eq 0 ]]; then
        touch metabat_dummy_${sample_id}.fa
    fi
    """
}

process CONCOCT {
    input:
    tuple val(sample_id), path(assembly), path(bam), path(bam_idx), path(bed)

    output:
    tuple val(sample_id), path('fasta_bins/*.fa'), emit: fasta, optional: true
    path 'contigs10k.fasta'
    path 'covtable.tsv'
    path 'clustering_gt1000.csv'
    path 'merged.csv'

    script:
    """
    cp "$bed" copy.bed # This is needed so Nextflow won't modify the input via the symlink

    cut_up_fasta.py "$assembly" -c 100000 --merge_last -b copy.bed > contigs10k.fasta
    
    concoct_coverage_table.py copy.bed "$bam" > covtable.tsv
    
    concoct --composition_file contigs10k.fasta --coverage_file covtable.tsv -s 100
    
    merge_cutup_clustering.py clustering_gt1000.csv > merged.csv

    mkdir -p fasta_bins
    extract_fasta_bins.py "$assembly" merged.csv --output_path fasta_bins
    """
}

process DREP {
    publishDir "${params.outdir}/drep", mode: "copy", pattern: "dereplicated_genomes"
    
    input:
    tuple val(sample_id), path(concoct_fa), path(maxbin_fa), path(metabat_fa)

    output:
    path 'drep_out/data_tables'
    path 'drep_out/dereplicated_genomes', emit: derepped
    path 'drep_out/data/checkM/checkM_outdir/results.tsv', emit: checkM_result

    script:
    """
    dRep dereplicate drep_out -g *.f*a
    """
}

process MAP2ASSEMBLY {
    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.bam"),
          path("${sample_id}.bam.bai"),
          path("${sample_id}.bed")

    script:
    """
    bwa index -p "$sample_id" "$assembly"

    bwa mem -t $task.cpus -a "$sample_id" ${reads[0]} ${reads[1]} |
        samtools sort -o "$sample_id".bam -

    samtools index "$sample_id".bam

    bedtools bamtobed -i "$sample_id".bam > "$sample_id".bed
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
    publishDir "${params.outdir}/metaphlan", mode: "copy", pattern: "*_profile.txt"
    publishDir "${params.outdir}/metaphlan", mode: "copy", pattern: "*.biom"

    input:
    tuple val(sample_id), path(reads)
    path metaphlan_db

    output:
    tuple val(sample_id), path("*_profile.txt")   ,                emit: profile
    tuple val(sample_id), path("*.biom")          ,                emit: biom
    tuple val(sample_id), path('*.bowtie2out.txt'), optional:true, emit: bt2out
    path "*_profile.txt"                                         , emit: mqc

    script:
    """
    BT2_DB_INDEX=`find -L ${metaphlan_db} -name "*.rev.1.bt2*" | sed 's/\\.rev.1.bt2.*\$//' | sed 's/.*\\///'`

    metaphlan \\
        --nproc ${task.cpus} \\
        --input_type fastq \\
        ${reads[0]},${reads[1]} \\
        --bowtie2out ${sample_id}.bowtie2out.txt \\
        --bowtie2db ${metaphlan_db} \\
        --index \$BT2_DB_INDEX \\
        --biom ${sample_id}.biom \\
        --output_file ${sample_id}_profile.txt
    """
}

process METAPHLAN_MERGE {
    input:
    path(profiles)

    output:
    path '*txt'

    script:
    """
    merge_metaphlan_tables.py \\
        -o ${sample_id}.txt \\
        ${profiles}
    """
}

//TODO FINISH THIS
process HUMANN {
    publishDir "${params.outdir}/humann", mode: "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path ('*{.log,.tsv}')

    script:
    """
    humann \
        --input $reads \
        --output . \
        --threads $task.cpus
    """
}
