#!/usr/bin/env nextflow

// Process parameters
krakendb_libs = params.krakendb_libs ? params.krakendb_libs?.split(',') as List : null

// Report
log.info """
    M E T A G E N O M I C S - N F   P I P E L I N E
    ==============================================================
    Reads in FASTQ files                                : ${params.reads}
    Output directory                                    : ${params.outdir}
    Path to pre-existing Kraken DB (if any)             : ${params.krakendb}
    Dir with genomes for custom Kraken DB (if any)      : ${params.krakendb_add}
    Library/ies for custom Kraken DB (if any)           : ${params.krakendb_libs}
    Path to pre-existing Kraken library dir (if any)    : ${params.krakendb_libdir}
    Taxonomic level for Bracken                         : ${params.bracken_taxlevel}
    Read length for Bracken                             : ${params.bracken_readlen}
    Min. nr. of reads for Bracken                       : ${params.bracken_minreads}
    ==============================================================
    """
    .stripIndent(true)

// Import processes
include { FASTP } from './modules.nf'
include { KRAKENDB_DL_DB } from './modules.nf'
include { KRAKENDB_DL_TAX } from './modules.nf'
include { KRAKENDB_DL_LIB } from './modules.nf'
include { KRAKENDB_COMBINE_AND_ADD } from './modules.nf'
include { KRAKENDB_BUILD } from './modules.nf'
include { KRAKENDB_COMBINE_LIBS } from './modules.nf'
include { KRAKEN } from './modules.nf'
include { BRACKENDB_BUILD } from './modules.nf'
include { BRACKEN } from './modules.nf'
include { ASSEMBLY } from './modules.nf'
include { MAP2ASSEMBLY } from './modules.nf'
include { MULTIQC } from './modules.nf'
include { HOST_INDEX } from './modules.nf'
include { HOST_REMOVE } from './modules.nf'
include { CONCOCT } from './modules.nf'
include { METABAT2 } from './modules.nf'
include { MAXBIN2 } from './modules.nf'
include { DREP } from './modules.nf'

// Define the workflow
workflow  {
    // Create channels from input files
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    host_asm_ch = Channel.fromPath(params.host_asm).first()
    krakendb_ch = params.krakendb
        ? Channel.fromPath(params.krakendb, checkIfExists: true).first()
        : null
    brackendb_ch = params.brackendb
        ? Channel.fromPath(params.brackendb, checkIfExists: true).first()
        : null
    krakendb_add_ch = params.krakendb_add
        ? Channel.fromPath(params.krakendb_add, checkIfExists: true).first()
        : Channel.empty()
    krakendb_liblist_ch = krakendb_libs ? Channel.fromList(krakendb_libs) : null
    krakendb_lib_ch = params.krakendb_libdir
        ? Channel.fromPath(params.krakendb_libdir, checkIfExists: true)
        : null
    
    // Create an index for the host reference genome, or use a pre-existing one
    host_index_ch = params.host_index
        ? Channel.fromPath(params.host_index, checkIfExists: true)
        : (HOST_INDEX(host_asm_ch))

    // Read preprocessing and QC
    fastp_ch = FASTP(reads_ch)
    
    // Kraken and Bracken DB-building
    if (!krakendb_ch) {
        krakendb_tax_ch = KRAKENDB_DL_TAX()
        // If no library-dir was provided, download libraries:
        if (!krakendb_lib_ch) {
            krakendb_lib_ch = KRAKENDB_DL_LIB(krakendb_liblist_ch).collect()
            krakendb_lib_ch = KRAKENDB_COMBINE_LIBS(krakendb_lib_ch)
        }
        // Combine taxonomy, libraries, and optionally custom-addition genomes:
        krakendb_unbuilt_ch = KRAKENDB_COMBINE_AND_ADD(
            krakendb_tax_ch,
            krakendb_lib_ch,
            krakendb_add_ch.ifEmpty(file('no_add'))
            )
        krakendb_ch = KRAKENDB_BUILD(krakendb_unbuilt_ch).first()
    }     
    
    if (!bracken_db_ch) {
        bracken_db_ch = BRACKENDB_BUILD(krakendb_ch, params.bracken_readlen)
    }

    // Run Kraken and Bracken
    kraken_ch = KRAKEN(fastp_ch.fastq, krakendb_ch)
    bracken_ch = BRACKEN(kraken_ch.report, bracken_db_ch,
                        params.bracken_taxlevel, params.bracken_minreads,
                        params.bracken_readlen)

    // Host removal
    host_remove_ch = HOST_REMOVE(host_index_ch, fastp_ch.fastq)
    //TODO Kraken after host removal? Or use kraken_extract to do host removal? 

    // Assembly
    asm_ch = ASSEMBLY(host_remove_ch)
    asm_map_ch = MAP2ASSEMBLY(fastp_ch.fastq, asm_ch.fasta)
    asm_and_map_ch = asm_ch.fasta.join(asm_map_ch)
    asm_and_reads_ch = asm_ch.fasta.join(fastp_ch.fastq)

    // binners
    concoct_ch = CONCOCT(asm_and_map_ch)
    maxbin_ch = MAXBIN2(asm_and_reads_ch)
    metabat_ch = METABAT2(asm_and_map_ch)
    bins_ch = concoct_ch.fasta.join(maxbin_ch.fasta).join(metabat_ch.fasta)

    // deprelicate and checkm
    drep_ch = DREP(bins_ch)

    // Assembly QC
    // checkm, Busco, etc

    // MultiQC
    mqc_in_ch = kraken_ch.report_path.mix(fastp_ch.report).collect()
    MULTIQC(mqc_in_ch)
}

// Report
workflow.onComplete {
    if (workflow.success) {
        log.info ("\nThe pipeline has finished successfully! Final outputs are in the $params.outdir dir.")
    } else {
        log.info ("\nThe pipeline encountered an error and did not finish successfully")
    }
}
