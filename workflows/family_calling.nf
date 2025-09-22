/*
 * Family Calling Workflow using GLnexus
 * Merges individual gVCF files into family VCF files
 */

// Include modules
include { GLNEXUS_FAMILY } from '../modules/glnexus_family'
include { NORMALIZE_VCF } from '../modules/normalize_vcf'

workflow FAMILY_CALLING {
    
    take:
    family_gvcfs    // channel: [fid, [barcodes], [gvcfs], [tbis]]
    
    main:
    
    // Run GLnexus for each family
    GLNEXUS_FAMILY(family_gvcfs)
    
    // Normalize the GLnexus output
    NORMALIZE_VCF(GLNEXUS_FAMILY.out.family_vcf)
    
    emit:
    family_vcfs = GLNEXUS_FAMILY.out.family_vcf
    normalized_vcfs = NORMALIZE_VCF.out.normalized_vcf
}