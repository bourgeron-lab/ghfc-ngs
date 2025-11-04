/*
 * GHFC WGS Alignment Workflow
 * Handles BWA-MEM2 alignment from FASTQ files, realignment from CRAM files,
 * and bedgraph generation from final CRAM files
 */

// Include modules
include { A_BWA_MEM2_ALIGN } from '../modules/a_bwa_mem2'
include { A_BAZAM_BWA_MEM2_REALIGN } from '../modules/a_bazam_bwa_mem2'
include { A_BAZAM_BWA_MEM2_REALIGN as A_BAZAM_BWA_MEM2_REALIGN_37 } from '../modules/a_bazam_bwa_mem2'
include { A_BAZAM_BWA_MEM2_REALIGN as A_BAZAM_BWA_MEM2_REALIGN_38 } from '../modules/a_bazam_bwa_mem2'
include { A_MERGE_UNITS } from '../modules/a_merge_units'
include { A_INDEX_CRAM } from '../modules/a_index_cram'
include { A_MOSDEPTH } from '../modules/a_mosdepth'
include { A_TABIX_INDEX } from '../modules/a_tabix_index'

workflow ALIGNMENT {
    
    take:
    fastq_files     // channel: [barcode, unit, fastq1, fastq2, project, flowcell, dual, lane]
    cram_37_files   // channel: [barcode, cram, crai]
    cram_38_files   // channel: [barcode, cram, crai]
    
    main:
    
    // BWA-MEM2 alignment from FASTQ files
    A_BWA_MEM2_ALIGN(
        fastq_files,
        params.ref,
        params.bwa_mem2,
        params.samblaster,
        params.sambamba,
        params.samtools,
        params.scratch
    )
    
    // Index intermediate CRAM files
    A_INDEX_CRAM(A_BWA_MEM2_ALIGN.out)
    
    // Group units by barcode for merging
    units_by_barcode = A_BWA_MEM2_ALIGN.out
        .map { barcode, unit, cram -> [barcode, cram] }
        .groupTuple()
        .join(
            A_INDEX_CRAM.out
                .map { barcode, unit, crai -> [barcode, crai] }
                .groupTuple()
        )
    
    // Merge units
    A_MERGE_UNITS(
        units_by_barcode,      // tuple val(barcode), path(crams), path(crais)
        params.samtools,  // val samtools_path
        params.data,      // val data_dir
        params.ref_name   // val ref_name
    )
    
    // Bazam realignment from GRCh37 CRAM files
    cram_37_with_oldref = cram_37_files.map { barcode, cram, crai -> 
        [barcode, cram, crai, params.old_ref_37] 
    }
    
    A_BAZAM_BWA_MEM2_REALIGN_37(
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
    
    A_BAZAM_BWA_MEM2_REALIGN_38(
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
    all_crams = A_MERGE_UNITS.out.cram
        .mix(A_BAZAM_BWA_MEM2_REALIGN_37.out.cram)
        .mix(A_BAZAM_BWA_MEM2_REALIGN_38.out.cram)

    // Generate bedgraphs from final CRAM files
    A_MOSDEPTH(
        all_crams,
        params.ref,
        params.bin,
        params.data
    )
    
    // Index bedgraph files with tabix
    A_TABIX_INDEX(
        A_MOSDEPTH.out.bedgraph,
        params.data
    )

    emit:
    crams = all_crams    // channel: [barcode, cram_file, crai_file]
    bedgraphs = A_TABIX_INDEX.out.indexed_bedgraph    // channel: [barcode, bedgraph, tbi]
}
