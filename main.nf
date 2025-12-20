#!/usr/bin/env nextflow

// Import processes
include { FASTQC } from './modules.nf'
include { FASTP } from './modules.nf'
include { READCOUNT as READCOUNT_PRE } from './modules.nf'
include { READCOUNT as READCOUNT_POST } from './modules.nf'
include { REPORT_READLOSS } from './modules.nf'
include { KRAKENDB_DL_DB } from './modules.nf'
include { KRAKENDB_DL_TAX } from './modules.nf'
include { KRAKENDB_DL_LIB } from './modules.nf'
include { KRAKENDB_COMBINE_AND_ADD } from './modules.nf'
include { KRAKENDB_BUILD } from './modules.nf'
include { KRAKENDB_COMBINE_LIBS } from './modules.nf'
include { BRACKENDB_BUILD } from './modules.nf'
include { KRAKEN as KRAKEN_HOSTRM } from './modules.nf'
include { KRAKEN as KRAKEN_CLASSIFY_READS_HIGH } from './modules.nf'
include { KRAKEN as KRAKEN_CLASSIFY_READS_LOW } from './modules.nf'
include { KRAKEN as KRAKEN_CLASSIFY_ASM_LOW } from './modules.nf'
include { KRAKEN as KRAKEN_CLASSIFY_ASM_HIGH } from './modules.nf'
include { KRAKEN_EXTRACT } from './modules.nf'
include { KRONA as KRONA_ON_HOSTRM } from './modules.nf'
include { KRONA as KRONA_ON_CLASSIFY_HIGH } from './modules.nf'
include { KRONA as KRONA_ON_CLASSIFY_LOW } from './modules.nf'
include { KRONA as KRONA_ON_BRACKEN_HIGH } from './modules.nf'
include { KRONA as KRONA_ON_BRACKEN_LOW } from './modules.nf'
include { KRONA_TAX } from './modules.nf'
include { BRACKEN as BRACKEN_HIGH } from './modules.nf'
include { BRACKEN as BRACKEN_LOW } from './modules.nf'
include { BIOM as KRAKEN_BIOM_HIGH } from './modules.nf'
include { BIOM as KRAKEN_BIOM_LOW } from './modules.nf'
include { BIOM as BRACKEN_BIOM_HIGH } from './modules.nf'
include { BIOM as BRACKEN_BIOM_LOW } from './modules.nf'
include { ASSEMBLY } from './modules.nf'
include { MAP2ASSEMBLY } from './modules.nf'
include { MULTIQC as MULTIQC_QC } from './modules.nf'
include { MULTIQC as MULTIQC_HOSTREMOVE } from './modules.nf'
include { MULTIQC as MULTIQC_READCLASSIF } from './modules.nf'
include { HOST_INDEX } from './modules.nf'
include { HOST_REMOVE_ALIGN } from './modules.nf'
include { CONCOCT } from './modules.nf'
include { METABAT2 } from './modules.nf'
include { MAXBIN2 } from './modules.nf'
include { DREP } from './modules.nf'
include { METAPHLAN_DB } from './modules.nf'
include { METAPHLAN } from './modules.nf'
include { METAPHLAN_MERGE } from './modules.nf'
include { SOURMASH_DB } from './modules.nf'
include { SOURMASH } from './modules.nf'

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
    conf_host = params.kraken_host_confidence
    minhit_host = params.kraken_host_minhitgroups

    minhit_classif_high = params.kraken_classif_minhitgroups_high
    conf_classif_high = params.kraken_classif_confidence_high
    minhit_classif_low = params.kraken_classif_minhitgroups_low
    conf_classif_low = params.kraken_classif_confidence_low

    // Report
    log.info """
    M E T A G E N O M I C S - N F   P I P E L I N E
    ==============================================================
    Reads in FASTQ files                   : ${params.reads}
    Output directory                       : ${params.outdir}
    Host removal method                    : ${params.host_removal_method}
    -----
    Host assembly genome (if any)          : ${params.host_asm}
    Kraken host-removal DB (if any)        : ${params.krakendb_host}
    Kraken assignment DB (if any)          : ${params.krakendb_classif}
    -----
    Skip assembly step?                    : ${params.skip_assembly}
    Skip Kraken classification?            : ${skip_kraken}
    Skip Bracken abundance estimation?     : ${skip_bracken}
    Skip DREP MAG dereplication?           : ${params.skip_drep}
    ==============================================================
    """.stripIndent(true)

    // =========================================================================
    //                   CREATE channelS FROM INPUT FILES
    // =========================================================================
    asm_ch = channel.empty()
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
    host_asm_ch = params.host_asm && ( params.host_removal_method == 'align' || params.host_removal_method == 'both' )
        ? channel.fromPath(params.host_asm).first()
        : channel.empty()
    metaphlandb_ch = params.metaphlandb
        ? channel.fromPath(params.metaphlandb, checkIfExists: true).first()
        : channel.empty()
    krakendb_host_ch = params.krakendb_host && ( params.host_removal_method == 'kraken' || params.host_removal_method == 'both' )
        ? channel.fromPath(params.krakendb_host, checkIfExists: true).first()
        : channel.empty()
    krakendb_classif_ch = params.krakendb_classif && !skip_kraken
        ? channel.fromPath(params.krakendb_classif, checkIfExists: true).first()
        : channel.empty()
    krakendb_host_add_ch = params.krakendb_host_add
        ? channel.fromPath(params.krakendb_host_add, checkIfExists: true).first()
        : channel.empty()
    krakendb_host_liblist_ch = krakendb_host_libs
        ? channel.fromList(krakendb_host_libs)
        : null
    krakendb_classif_add_ch = params.krakendb_classif_add
        ? channel.fromPath(params.krakendb_classif_add, checkIfExists: true).first()
        : channel.empty()
    krakendb_classif_liblist_ch = krakendb_classif_libs
        ? channel.fromList(krakendb_classif_libs)
        : null
    brackendb_ch = params.brackendb && !skip_bracken
        ? channel.fromPath(params.brackendb, checkIfExists: true).first()
        : channel.empty()
    kraken_host_ch = channel.empty()
    host_aln_ch = channel.empty()

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
            krakendb_host_add_ch.ifEmpty(file('no_add'))
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
            krakendb_classif_add_ch.ifEmpty(file('no_add'))
            )
        krakendb_classif_ch = KRAKENDB_BUILD(krakendb_unbuilt_ch).first()
    }
    
    // =========================================================================
    //                  BASIC READ QC AND PREPROCESSING
    // =========================================================================
    // FastQC
    fastqc_ch = FASTQC(reads_ch)
    
    // Fastp
    fastp_ch = FASTP(reads_ch)
    reads_ch = fastp_ch.fastq

    // Pre-host removal read count - the input should be a simple list of FASTQ files:
    fastq_list_ch = reads_ch.map {_id, reads -> reads[0]}.flatten().collect()
    readcount_pre_ch = READCOUNT_PRE(fastq_list_ch, 'pre-host-remove')

    // =========================================================================
    //                  HOST READ REMOVAL
    // =========================================================================
    // Alignmnent-based:
    if (params.host_removal_method == 'align' || params.host_removal_method == 'both') {
        
        // Create an index for the host reference genome, or use a pre-existing one
        host_index_ch = params.host_index
            ? channel.fromPath(params.host_index, checkIfExists: true).collect()
            : HOST_INDEX(host_asm_ch)
        
        // Align the reads
        host_aln_ch = HOST_REMOVE_ALIGN(host_index_ch, reads_ch)
        reads_ch = host_aln_ch.fastq
    }

    // Kraken-based:
    if (params.host_removal_method == 'kraken' || params.host_removal_method == 'both') {
        
        // Host read removal with Kraken
        kraken_host_ch = KRAKEN_HOSTRM(
            reads_ch, krakendb_host_ch, conf_host, minhit_host, 'hostremove'
        )

        // Extract non-assigned reads
        extract_input_ch = kraken_host_ch.k_extract.join(reads_ch)
        reads_ch = KRAKEN_EXTRACT(extract_input_ch, params.kraken_tax_remove).fastq
        
        // Krona
        KRONA_ON_HOSTRM(kraken_host_ch.main_out, krona_tax_sh)
    }

    // Post-host removal read count - the input should be a simple list of FASTQ files:
    fastq_list_ch = reads_ch.map {_id, reads -> reads[0]}.flatten().collect()
    readcount_post_ch = READCOUNT_POST(fastq_list_ch, 'post-host-remove')
    REPORT_READLOSS(readcount_pre_ch, readcount_post_ch)

    // =========================================================================
    //                          READ CLASSIFICATION
    // =========================================================================
    // Kraken
    kraken_classif_high_ch = KRAKEN_CLASSIFY_READS_HIGH(
        reads_ch, krakendb_classif_ch, conf_classif_high, minhit_classif_high, 'classif_high'
        )
    kraken_classif_low_ch = KRAKEN_CLASSIFY_READS_LOW(
        reads_ch, krakendb_classif_ch, conf_classif_low, minhit_classif_low, 'classif_low'
        )
    KRONA_ON_CLASSIFY_HIGH(kraken_classif_high_ch.main_out, krona_tax_sh)
    KRONA_ON_CLASSIFY_LOW(kraken_classif_low_ch.main_out, krona_tax_sh)
    
    KRAKEN_BIOM_HIGH(kraken_classif_high_ch.mqc.collect(), 'kraken_high')
    KRAKEN_BIOM_LOW(kraken_classif_low_ch.mqc.collect(), 'kraken_low')

    // Bracken
    if (!brackendb_ch && !skip_bracken) {
        brackendb_ch = BRACKENDB_BUILD(krakendb_classif_ch, params.bracken_readlen)
    }

    bracken_high_ch = BRACKEN_HIGH(
        kraken_classif_high_ch.report, brackendb_ch,
        params.bracken_taxlevel, params.bracken_minreads, params.bracken_readlen, 'high'
    )
    // KRONA_ON_BRACKEN(bracken_ch.main_out, krona_tax_sh) // THIS WILL NOT WORK, NEED 'MAIN' KRAKEN-STYLE OUTPUT
    BRACKEN_BIOM_HIGH(bracken_high_ch.report.collect(), 'bracken')

    bracken_low_ch = BRACKEN_LOW(
        kraken_classif_low_ch.report, brackendb_ch,
        params.bracken_taxlevel, params.bracken_minreads, params.bracken_readlen, 'low'
    )
    // KRONA_ON_BRACKEN(bracken_ch.main_out, krona_tax_sh) // THIS WILL NOT WORK, NEED 'MAIN' KRAKEN-STYLE OUTPUT
    BRACKEN_BIOM_LOW(bracken_low_ch.report.collect(), 'bracken')

    // MetaPhlAn
    metaphlan_ch = METAPHLAN(reads_ch, metaphlandb_ch)
    METAPHLAN_MERGE(metaphlan_ch.mqc.collect())

    // =========================================================================
    //                              MAG ASSEMBLY
    // =========================================================================
    // Assembly
    if (params.skip_assembly == false) asm_ch = ASSEMBLY(reads_ch).fasta
    asm_and_reads_ch = asm_ch.join(reads_ch)
    asm_map_ch = MAP2ASSEMBLY(asm_and_reads_ch)
    asm_and_map_ch = asm_ch.join(asm_map_ch)

    // Binning
    // maxbin_ch = MAXBIN2(asm_and_reads_ch) -- DISABLED FOR NOW DUE TO ERRORS
    metabat_ch = METABAT2(asm_and_map_ch)
    concoct_ch = CONCOCT(asm_and_map_ch)
    // bins_ch = concoct_ch.fasta.join(maxbin_ch.fasta).join(metabat_ch.fasta)
    bins_ch = concoct_ch.fasta.join(metabat_ch.fasta)
    drep_ch = DREP(bins_ch)

    // =========================================================================
    //                              MAG CLASSIFICATION
    // =========================================================================
    // Sourmash
    //sourmash_db_ch = SOURMASH_DB()
    //SOURMASH(bins_ch, sourmash_db_ch)

    // Kraken
    //kraken_classif_high_ch = KRAKEN_CLASSIFY_ASM_HIGH(
    //    bins_ch, krakendb_classif_ch, conf_classif_high, minhit_classif_high, 'classif_high'
    //)
    //kraken_classif_low_ch = KRAKEN_CLASSIFY_ASM_LOW(
    //    bins_ch, krakendb_classif_ch, conf_classif_low, minhit_classif_low, 'classif_low'
    //)

    // =========================================================================
    //                              MULTIQC
    // =========================================================================
    // MultiQC
    mqc_qc_ch = fastqc_ch.zip.mix(fastp_ch.report).collect()
    MULTIQC_QC(mqc_qc_ch, 'read-qc')

    mqc_host_ch = host_aln_ch.logs.ifEmpty([])
        .mix(kraken_host_ch.mqc.ifEmpty([]))
        .collect()
    MULTIQC_HOSTREMOVE(mqc_host_ch, 'host-remove')

    mqc_readclassif_ch = kraken_classif_high_ch.mqc.ifEmpty([])
        .mix(kraken_classif_low_ch.mqc.ifEmpty([]))
        .mix(bracken_high_ch.report.ifEmpty([]))
        .mix(bracken_low_ch.report.ifEmpty([]))
        .mix(metaphlan_ch.mqc.ifEmpty([]))
        .collect()
    MULTIQC_READCLASSIF(mqc_readclassif_ch, 'read-classif')

    // =========================================================================
    //                              MULTIQC
    // =========================================================================
    workflow.onComplete = {
        log.info "\n===================\nPipeline completed at: $workflow.complete"
        log.info "\nThe pipeline ${ workflow.success ? 'completed successfully!' : 'failed!' }"
        log.info "=================="
    }

}
