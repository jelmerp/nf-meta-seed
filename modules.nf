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

    script:
    """
    mkdir -p lib_dir
    cp -rv $lib_dirs lib_dir/
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

    script:
    """
    mkdir -p kraken_db_unbuilt/library
    cp -rv ${tax_dir}/taxonomy kraken_db_unbuilt
    cp -rv ${lib_dir}/* kraken_db_unbuilt/library/

    if [[ -d ${genomes_to_add_dir} ]]; then
        for fasta in ${genomes_to_add_dir}/*; do
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

    script:
    """
    mkdir kraken_db
    cp -rv $kraken_db_unbuilt/* kraken_db/
    
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

    script:
    """
    cp -rv $kraken_db bracken_db

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
    path "${sample_id}.fasta", emit: assembly
    path "${sample_id}_scaffolds.fasta"
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
    
    mv outdir/contigs.fasta ${sample_id}.fasta
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
