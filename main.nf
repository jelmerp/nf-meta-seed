#!/usr/bin/env nextflow

// Import processes
include { FASTQC } from './modules.nf'
include { FASTP } from './modules.nf'
include { KRAKENDB_DL_DB } from './modules.nf'
include { KRAKENDB_DL_TAX } from './modules.nf'
include { KRAKENDB_DL_LIB } from './modules.nf'
include { KRAKENDB_COMBINE_AND_ADD } from './modules.nf'
include { KRAKENDB_BUILD } from './modules.nf'
include { KRAKENDB_COMBINE_LIBS } from './modules.nf'
include { BRACKENDB_BUILD } from './modules.nf'
include { KRAKEN as KRAKEN_HOSTRM } from './modules.nf'
include { KRAKEN as KRAKEN_CLASSIFY } from './modules.nf'
include { KRAKEN_EXTRACT } from './modules.nf'
include { KRONA as KRONA_ON_HOSTRM } from './modules.nf'
include { KRONA as KRONA_ON_CLASSIFY } from './modules.nf'
include { KRONA_TAX } from './modules.nf'
include { BRACKEN } from './modules.nf'
include { ASSEMBLY } from './modules.nf'
include { MAP2ASSEMBLY } from './modules.nf'
include { MULTIQC } from './modules.nf'
include { HOST_INDEX } from './modules.nf'
include { HOST_REMOVE_ALIGN } from './modules.nf'
include { CONCOCT } from './modules.nf'
include { METABAT2 } from './modules.nf'
include { MAXBIN2 } from './modules.nf'
include { DREP } from './modules.nf'
include { METAPHLAN_DB } from './modules.nf'
include { METAPHLAN } from './modules.nf'
include { METAPHLAN_MERGE } from './modules.nf'

// Define the workflow
workflow  {
    // Process parameters
    krakendb_host_libs = params.krakendb_host_libs
        ? params.krakendb_libs?.split(',') as List
        : null
    krakendb_classif_libs = params.krakendb_classif_libs
        ? params.krakendb_libs?.split(',') as List
        : null
    skip_kraken = params.skip_kraken
    skip_bracken = params.skip_bracken
    skip_bracken = skip_kraken ? true : skip_bracken
    conf_classif = params.kraken_classif_confidence
    conf_host = params.kraken_host_confidence
    minhit_classif = params.kraken_classif_minhitgroups
    minhit_host = params.kraken_host_minhitgroups

    // Report
    log.info """
    M E T A G E N O M I C S - N F   P I P E L I N E
    ==============================================================
    Reads in FASTQ files                   : ${params.reads}
    Output directory                       : ${params.outdir}
    Host removal method                    : ${params.host_removal_method}
    Host assembly genome (if any)          : ${params.host_asm}
    ==============================================================
    """.stripIndent(true)

    // =========================================================================
    //                   CREATE CHANNELS FROM INPUT FILES
    // =========================================================================
    asm_ch = Channel.empty()
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    host_asm_ch = params.host_asm && params.host_removal_method == 'align'
        ? Channel.fromPath(params.host_asm).first()
        : Channel.empty()
    metaphlandb_ch = params.metaphlandb
        ? Channel.fromPath(params.metaphlandb, checkIfExists: true).first()
        : Channel.empty()
    krakendb_host_ch = params.krakendb_host && params.host_removal_method == 'kraken'
        ? Channel.fromPath(params.krakendb_host, checkIfExists: true).first()
        : Channel.empty()
    krakendb_classif_ch = params.krakendb_classif && !skip_kraken
        ? Channel.fromPath(params.krakendb_classif, checkIfExists: true).first()
        : Channel.empty()
    krakendb_host_genomes_ch = params.krakendb_host_genomes
        ? Channel.fromPath(params.krakendb_host_genomes, checkIfExists: true).first()
        : Channel.empty()
    krakendb_host_liblist_ch = krakendb_host_libs
        ? Channel.fromList(krakendb_host_libs)
        : null
    krakendb_classif_genomes_ch = params.krakendb_classif_genomes
        ? Channel.fromPath(params.krakendb_classif_genomes, checkIfExists: true).first()
        : Channel.empty()
    krakendb_classif_liblist_ch = krakendb_classif_libs
        ? Channel.fromList(krakendb_classif_libs)
        : null
    brackendb_ch = params.brackendb && !skip_bracken
        ? Channel.fromPath(params.brackendb, checkIfExists: true).first()
        : Channel.empty()
    kraken_host_mqc_ch = Channel.empty()

    // =========================================================================
    //                          DATABASE BUILDING
    // =========================================================================
    // Krona
    krona_tax_sh = KRONA_TAX()
    
    // Metaphlan
    if (!metaphlandb_ch) metaphlandb_ch = METAPHLAN_DB()
    
    // Kraken db for host-removal
    if (!krakendb_host_ch && params.host_removal_method == 'kraken') {
        
        krakendb_host_tax_ch = KRAKENDB_DL_TAX()
        
        // If no library-dir was provided, download libraries:
        krakendb_host_lib_ch = KRAKENDB_DL_LIB(krakendb_host_liblist_ch).collect()
        krakendb_host_lib_ch = KRAKENDB_COMBINE_LIBS(krakendb_host_lib_ch)
        
        // Combine taxonomy, libraries, and optionally custom-addition genomes:
        krakendb_host_unbuilt_ch = KRAKENDB_COMBINE_AND_ADD(
            krakendb_host_tax_ch,
            krakendb_host_lib_ch,
            krakendb_host_genomes_ch.ifEmpty(file('no_add'))
            )
        krakendb_host_ch = KRAKENDB_BUILD(krakendb_host_unbuilt_ch).first()
    }

    // Classification Kraken db
    if (!krakendb_classif_ch && !skip_kraken) {
        krakendb_classif_tax_ch = KRAKENDB_DL_TAX()
        
        // If no library-dir was provided, download libraries:
        krakendb_classif_lib_ch = KRAKENDB_DL_LIB(krakendb_classif_liblist_ch).collect()
        krakendb_classif_lib_ch = KRAKENDB_COMBINE_LIBS(krakendb_classif_lib_ch)
        
        // Combine taxonomy, libraries, and optionally custom-addition genomes:
        krakendb_unbuilt_ch = KRAKENDB_COMBINE_AND_ADD(
            krakendb_classif_tax_ch,
            krakendb_classif_lib_ch,
            krakendb_classif_genomes_ch.ifEmpty(file('no_add'))
            )
        krakendb_classif_ch = KRAKENDB_BUILD(krakendb_unbuilt_ch).first()
    }
    
    // =========================================================================
    //                  READ QC AND PREPROCESSING
    // =========================================================================
    // FastQC
    fastqc_ch = FASTQC(reads_ch)
    
    // Fastp
    fastp_ch = FASTP(reads_ch)
    reads_ch = fastp_ch.fastq

    // Host read removal
    if (params.host_removal_method == 'kraken') {
        // Host read removal with Kraken
        kraken_host_ch = KRAKEN_HOSTRM(
            reads_ch, krakendb_host_ch, conf_host, minhit_host, 'hostremove'
            )
        kraken_host_mqc_ch = kraken_host_ch.mqc 
        extract_input_ch = kraken_host_ch.output.join(reads_ch)
        reads_ch = KRAKEN_EXTRACT(extract_input_ch, params.kraken_tax_remove).fq
        //TODO report how many reads were removed
        KRONA_ON_HOSTRM(kraken_host_ch.output, krona_tax_sh)
    
    } else if (params.host_removal_method == 'align') {
        // Create an index for the host reference genome, or use a pre-existing one
        host_index_ch = params.host_index
            ? Channel.fromPath(params.host_index, checkIfExists: true)
            : HOST_INDEX(host_asm_ch)
        reads_ch = HOST_REMOVE_ALIGN(host_index_ch, reads_ch).fastq
    }
    
    // =========================================================================
    //                          READ CLASSIFICATION
    // =========================================================================
    // Kraken
    kraken_classif_ch = KRAKEN_CLASSIFY(
        reads_ch, krakendb_classif_ch, conf_classif, minhit_classif, 'classify'
        )
    KRONA_ON_CLASSIFY(kraken_classif_ch.output, krona_tax_sh)
    
    // Bracken
    if (!brackendb_ch && !skip_bracken) {
        brackendb_ch = BRACKENDB_BUILD(krakendb_classif_ch, params.bracken_readlen)
    }
    bracken_ch = BRACKEN(
        kraken_classif_ch.output, brackendb_ch,
        params.bracken_taxlevel, params.bracken_minreads, params.bracken_readlen
        )
    //TODO - Krona on Bracken output?
    //TODO - Use https://github.com/jenniferlu717/KrakenTools?tab=readme-ov-file#kreport2kronapy

    // MetaPhlAn
    metaphlan_ch = METAPHLAN(reads_ch, metaphlandb_ch)
    //TODO - METAPHLAN_MERGE()
    //TODO - Strainphlan - https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1
    //TODO - Visualization with Graphphlan? https://github.com/biobakery/graphlan/wiki
    //TODO - PhyloPhlAn?

    // =========================================================================
    //                              MAG ASSEMBLY
    // =========================================================================
    // Assembly
    if (params.skip_assembly == false) asm_ch = ASSEMBLY(reads_ch).fasta
    asm_and_reads_ch = asm_ch.join(reads_ch)
    asm_map_ch = MAP2ASSEMBLY(asm_and_reads_ch)
    asm_and_map_ch = asm_ch.join(asm_map_ch)

    // Binning
    maxbin_ch = MAXBIN2(asm_and_reads_ch)
    metabat_ch = METABAT2(asm_and_map_ch)
    concoct_ch = CONCOCT(asm_and_map_ch)
    bins_ch = concoct_ch.fasta.join(maxbin_ch.fasta).join(metabat_ch.fasta)
    drep_ch = DREP(bins_ch)

    // TODO - Assembly QC - Busco, etc
    // TODO - Classification of the MAGs
    // TODO - Abundance estimation
    // TODO - Functional analysis

    // =========================================================================
    //                              MULTIQC
    // =========================================================================
    // MultiQC
    mqc_in_ch = fastqc_ch.zip
        .mix(fastp_ch.report)
        .mix(kraken_host_mqc_ch.ifEmpty([]))
        .mix(kraken_classif_ch.mqc.ifEmpty([]))
        .mix(bracken_ch.ifEmpty([]))
        .mix(metaphlan_ch.mqc.ifEmpty([]))
        .flatten()
        .collect()
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
