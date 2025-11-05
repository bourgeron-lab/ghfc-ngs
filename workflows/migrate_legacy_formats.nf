#!/usr/bin/env nextflow

/*
========================================================================================
    GHFC WGS Legacy Format Migration Workflow
========================================================================================
    One-time migration workflow to convert legacy VCF.gz files to BCF format
    and clean up obsolete intermediate files.
    
    This should be run BEFORE the main pipeline if you have legacy data.
========================================================================================
*/

nextflow.enable.dsl = 2

// Include migration modules
include { CONVERT_NORM_VCF_TO_BCF } from '../modules/migration/convert_norm_vcf_to_bcf'
include { CONVERT_COMMON_VCF_TO_BCF } from '../modules/migration/convert_common_vcf_to_bcf'
include { CLEANUP_LEGACY_FILES } from '../modules/migration/cleanup_legacy_files'

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.data) {
    exit 1, "ERROR: --data parameter is required"
}

log.info """
========================================================================================
                    GHFC WGS Legacy Format Migration
========================================================================================
Data directory   : ${params.data}
========================================================================================
This workflow will:
1. Convert {FID}.norm.vcf.gz → {FID}.norm.bcf
2. Convert {FID}.common.vcf.gz → {FID}.common.bcf
3. Remove obsolete intermediate files:
   - {FID}.vcf.gz (GLnexus output)
   - {FID}.gnomad.vcf.gz / {FID}.gnomad.bcf (gnomAD annotation)
========================================================================================
"""

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow MIGRATE_LEGACY_FORMATS {
    
    main:
    
    // Discover families from directory structure
    def families_dir = new File("${params.data}/families")
    if (!families_dir.exists()) {
        log.warn "No families directory found at: ${params.data}/families"
        log.warn "Nothing to migrate."
        return
    }
    
    def families = []
    families_dir.eachDir { family_dir ->
        families.add(family_dir.name)
    }
    
    if (families.isEmpty()) {
        log.warn "No family directories found in: ${params.data}/families"
        log.warn "Nothing to migrate."
        return
    }
    
    log.info "Found ${families.size()} families to check for migration"
    
    // Check for legacy norm.vcf.gz files and create channel
    def norm_files_to_convert = []
    families.each { fid ->
        def norm_vcf = new File("${params.data}/families/${fid}/vcfs/${fid}.norm.vcf.gz")
        def norm_tbi = new File("${params.data}/families/${fid}/vcfs/${fid}.norm.vcf.gz.tbi")
        def norm_bcf = new File("${params.data}/families/${fid}/vcfs/${fid}.norm.bcf")
        
        if (norm_vcf.exists() && norm_tbi.exists() && !norm_bcf.exists()) {
            log.info "  Found legacy norm.vcf.gz for family ${fid}"
            norm_files_to_convert.add([fid, norm_vcf, norm_tbi])
        }
    }
    
    // Check for legacy common.vcf.gz files and create channel
    def common_files_to_convert = []
    families.each { fid ->
        def common_vcf = new File("${params.data}/families/${fid}/vcfs/${fid}.common.vcf.gz")
        def common_tbi = new File("${params.data}/families/${fid}/vcfs/${fid}.common.vcf.gz.tbi")
        def common_bcf = new File("${params.data}/families/${fid}/vcfs/${fid}.common.bcf")
        
        if (common_vcf.exists() && common_tbi.exists() && !common_bcf.exists()) {
            log.info "  Found legacy common.vcf.gz for family ${fid}"
            common_files_to_convert.add([fid, common_vcf, common_tbi])
        }
    }
    
    log.info ""
    log.info "Migration plan:"
    log.info "  - ${norm_files_to_convert.size()} norm.vcf.gz files to convert"
    log.info "  - ${common_files_to_convert.size()} common.vcf.gz files to convert"
    log.info "  - ${families.size()} families to clean up"
    log.info ""
    
    // Convert norm.vcf.gz files to BCF
    if (!norm_files_to_convert.isEmpty()) {
        norm_channel = Channel.fromList(norm_files_to_convert)
            .map { fid, vcf, tbi -> [fid, file(vcf), file(tbi)] }
        
        CONVERT_NORM_VCF_TO_BCF(norm_channel)
    }
    
    // Convert common.vcf.gz files to BCF
    if (!common_files_to_convert.isEmpty()) {
        common_channel = Channel.fromList(common_files_to_convert)
            .map { fid, vcf, tbi -> [fid, file(vcf), file(tbi)] }
        
        CONVERT_COMMON_VCF_TO_BCF(common_channel)
    }
    
    // Clean up legacy files for all families
    cleanup_channel = Channel.fromList(families)
        .map { fid -> [fid, params.data, []] }
    
    CLEANUP_LEGACY_FILES(cleanup_channel)
    
    emit:
    norm_converted = norm_files_to_convert.isEmpty() ? Channel.empty() : CONVERT_NORM_VCF_TO_BCF.out.converted_bcf
    common_converted = common_files_to_convert.isEmpty() ? Channel.empty() : CONVERT_COMMON_VCF_TO_BCF.out.converted_bcf
    cleanup_status = CLEANUP_LEGACY_FILES.out.status
}

/*
========================================================================================
    COMPLETION HANDLER
========================================================================================
*/

workflow.onComplete {
    log.info ""
    log.info "========================================================================================"
    log.info "Migration workflow completed!"
    log.info "========================================================================================"
    log.info "Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration: ${workflow.duration}"
    log.info ""
    if (workflow.success) {
        log.info "✓ Legacy VCF.gz files have been converted to BCF format"
        log.info "✓ Obsolete intermediate files have been removed"
        log.info ""
        log.info "You can now run the main pipeline safely."
    } else {
        log.info "✗ Migration failed. Please check the error messages above."
    }
    log.info "========================================================================================"
}

/*
========================================================================================
    RUN WORKFLOW
========================================================================================
*/

workflow {
    MIGRATE_LEGACY_FORMATS()
}
