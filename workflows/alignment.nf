/*
 * GHFC WGS Alignment Workflow
 * Handles BWA-MEM2 alignment from FASTQ files and realignment from CRAM files
 */

// Include modules
include { BWA_MEM2_ALIGN } from '../modules/bwa_mem2'
include { BAZAM_BWA_MEM2_REALIGN } from '../modules/bazam_bwa_mem2'
include { MERGE_UNITS } from '../modules/merge_units'
include { INDEX_CRAM } from '../modules/index_cram'

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
    
    BAZAM_BWA_MEM2_REALIGN(
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
    
    BAZAM_BWA_MEM2_REALIGN(
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
    all_crams = MERGE_UNITS.out.cram.mix(BAZAM_BWA_MEM2_REALIGN.out.cram)

    emit:
    crams = all_crams    // channel: [barcode, cram_file, crai_file]
}
