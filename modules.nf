// QC and trim reads
process FASTP {
    publishDir "$params.outdir/fastp", mode: 'copy', pattern: '*html'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*fastq.gz'), emit: fastq
    path '*fastp_report*', emit: report

    script:
    """
    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -O ${sample_id}_R2_trimmed.fastq.gz \
        -h ${sample_id}_fastp_report.html \
        -j ${sample_id}_fastp_report.json \
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
    publishDir "$params.outdir/kraken", mode: 'copy', pattern: '*txt'

    input:
    tuple val(sample_id), path(reads)
    path(db)

    output:
    tuple val(sample_id), path('*fastq'), emit: classified_fq
    tuple val(sample_id), path('*report.txt'), emit: report
    tuple val(sample_id), path('*main.txt'), emit: main
    path '*report.txt', emit: report_path // For MultiQC, to avoid complications with tuple

    when:
    params.skip_kraken == false

    script:
    """
    kraken2 \
        --db ${db} \
        --report ${sample_id}_kraken-report.txt \
        --output ${sample_id}_kraken-main.txt \
        --classified-out ${sample_id}#.fastq \
        --gzip-compressed \
        --paired \
        --threads $task.cpus \
        ${reads[0]} \
        ${reads[1]}
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

// Run Bracken
process BRACKEN {
    publishDir "$params.outdir/bracken", mode: 'copy', pattern: '*txt'

    input:
    tuple val(sample_id), path(kraken_report)
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
        -i ${kraken_report} \
        -d ${bracken_db} \
        -o ${sample_id}_bracken-out.txt \
        -w ${sample_id}_bracken-report.txt \
        -r ${read_len} \
        -l ${tax_level} \
        -t ${min_reads}
    """
}

// Assemble reads
process ASSEMBLY {
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_scaffolds.fasta"), emit: fasta
    path "${sample_id}_contigs.fasta"
    path "${sample_id}_spades.log"
    
    script:
    def memory_gb = MemoryUnit.of("${task.memory}").toUnit('GB')
    """
    spades.py \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o outdir \
        --only-assembler \
        --meta \
        --threads $task.cpus \
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
    multiqc .
    """
}

process HOST_INDEX {
    publishDir "${params.outdir}/hostindex", mode: "copy"
    
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

process HOST_REMOVE {
    input:
    path(host_index_dir)
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*.fastq')

    shell:
    '''
    index_prefix=$(ls !{host_index_dir} | head -n1 | sed -E "s/.[0-9]+.bt2//")
    index_prefix_full=!{host_index_dir}/$index_prefix

    bowtie2 \
        -p !{task.cpus} \
        -x $index_prefix_full \
        -1 !{reads[0]} \
        -2 !{reads[1]} \
        --local \
        --un-conc \
        !{sample_id}_host_removed_reads \
        > !{sample_id}_mapped_unmapped.sam

    mv !{sample_id}_host_removed_reads.1 !{sample_id}_hostrm_R1.fastq
    mv !{sample_id}_host_removed_reads.2 !{sample_id}_hostrm_R2.fastq
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
    
    run_MaxBin.pl \
        -contig $assembly \
        -reads ${reads[0]} \
        -reads2 ${reads[1]} \
        -out $sample_id

    [[ \$? -ne 0 ]] && touch maxbin_dummy_${sample_id}.fasta
    """
}

process METABAT2 {
    input:
    tuple val(sample_id), path(assembly), path(bam), path(bam_idx), path(bed)

    output:
    tuple val(sample_id), path ('*.fasta'), emit: fasta, optional: true
    path 'depth.txt'

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt "$bam"

    metabat2 -i "$assembly" -a depth.txt -m 1500 --maxP 75 -s 100000 -o "$sample_id"

    n_files=`find . -type f -name "*fasta" | wc -l`
    [[ \$n_files -eq 0 ]] && touch metabat_dummy_${sample_id}.fasta
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
    cut_up_fasta.py "$assembly" -c 100000 --merge_last -b "$bed" > contigs10k.fasta
    concoct_coverage_table.py "$bed" "$bam" > covtable.tsv
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
    path 'data_tables'
    path 'dereplicated_genomes', emit: depreplicated_genomes
    path 'data/checkM/checkM_outdir/results.tsv', emit: checkM_result

    script:
    """
    dRep dereplicate drep_out -g *.f*a
    """
}

process MAP2ASSEMBLY {
    input:
    tuple val(sample_id), path(reads)
    tuple val(sample_id), path(assembly)

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
