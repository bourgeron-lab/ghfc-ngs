/*
 * GHFC WGS Alignment Workflow
 * Handles BWA-MEM2 alignment from FASTQ files, realignment from CRAM files,
 * and bedgraph generation from final CRAM files
 */

// Include modules
include { BWA_MEM2_ALIGN } from '../modules/alignment/bwa_mem2'
include { BAZAM_BWA_MEM2_REALIGN } from '../modules/alignment/bazam_bwa_mem2'
include { BAZAM_BWA_MEM2_REALIGN as BAZAM_BWA_MEM2_REALIGN_37 } from '../modules/alignment/bazam_bwa_mem2'
include { BAZAM_BWA_MEM2_REALIGN as BAZAM_BWA_MEM2_REALIGN_38 } from '../modules/alignment/bazam_bwa_mem2'
include { MERGE_UNITS } from '../modules/alignment/merge_units'
include { INDEX_CRAM } from '../modules/alignment/index_cram'
include { MOSDEPTH } from '../modules/alignment/mosdepth'
include { TABIX_INDEX } from '../modules/alignment/tabix_index'

workflow ALIGNMENT {
    
    take:
    fastq_files     // channel: [barcode, unit, fastq1, fastq2, project, flowcell, dual, lane]
    cram_37_files   // channel: [barcode, cram, crai]
    cram_38_files   // channel: [barcode, cram, crai]
    
    main:
    
    // BWA-MEM2 alignment from FASTQ files
    BWA_MEM2_ALIGN(
        fastq_files,
        params.ref,
        params.bwa_mem2,
        params.samblaster,
        params.sambamba,
        params.samtools,
        params.scratch
    )
    
    // Index intermediate CRAM files
    INDEX_CRAM(BWA_MEM2_ALIGN.out)
    
    // Group units by barcode for merging
    units_by_barcode = BWA_MEM2_ALIGN.out
        .map { barcode, unit, cram -> [barcode, cram] }
        .groupTuple()
        .join(
            INDEX_CRAM.out
                .map { barcode, unit, crai -> [barcode, crai] }
                .groupTuple()
        )
    
    // Merge units
    MERGE_UNITS(
        units_by_barcode,      // tuple val(barcode), path(crams), path(crais)
        params.samtools,  // val samtools_path
        params.data,      // val data_dir
        params.ref_name   // val ref_name
    )
    
    // Bazam realignment from GRCh37 CRAM files
    cram_37_with_oldref = cram_37_files.map { barcode, cram, crai -> 
        [barcode, cram, crai, params.old_ref_37] 
    }
    
    BAZAM_BWA_MEM2_REALIGN_37(
        cram_37_with_oldref,
        params.ref,
        params.bwa_mem2,
        params.bazam,
        params.samblaster,
        params.sambamba,
        params.samtools,
        params.scratch,
        params.data,
        params.ref_name
    )
    
    // Bazam realignment from GRCh38 CRAM files
    cram_38_with_oldref = cram_38_files.map { barcode, cram, crai -> 
        [barcode, cram, crai, params.old_ref_38] 
    }
    
    BAZAM_BWA_MEM2_REALIGN_38(
        cram_38_with_oldref,
        params.ref,
        params.bwa_mem2,
        params.bazam,
        params.samblaster,
        params.sambamba,
        params.samtools,
        params.scratch,
        params.data,
        params.ref_name
    )
    
    // Combine all final CRAM outputs
    all_crams = MERGE_UNITS.out.cram
        .mix(BAZAM_BWA_MEM2_REALIGN_37.out.cram)
        .mix(BAZAM_BWA_MEM2_REALIGN_38.out.cram)

    // Generate bedgraphs from final CRAM files
    MOSDEPTH(
        all_crams,
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
    crams = all_crams    // channel: [barcode, cram_file, crai_file]
    bedgraphs = TABIX_INDEX.out.indexed_bedgraph    // channel: [barcode, bedgraph, tbi]
}
