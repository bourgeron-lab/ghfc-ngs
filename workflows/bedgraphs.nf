/*
 * GHFC WGS Bedgraphs Workflow
 * Generates bedgraph files from CRAM files using mosdepth and indexes them with tabix
 */

// Include modules
include { MOSDEPTH } from '../modules/mosdepth'
include { TABIX_INDEX } from '../modules/tabix_index'

workflow BEDGRAPHS {
    
    take:
    cram_files     // channel: [barcode, cram, crai]
    
    main:
    
    // Run mosdepth to generate bedgraph files
    MOSDEPTH(
        cram_files,
        params.ref,
        params.bin,
        params.data
    )
    
    // Index bedgraph files with tabix
    TABIX_INDEX(
        MOSDEPTH.out.bedgraph,
        params.data
    )

    emit:
    bedgraphs = TABIX_INDEX.out.indexed_bedgraph    // channel: [barcode, bedgraph, tbi]
}
