/*
 * Family Calling Workflow using GLnexus
 * Merges individual gVCF files into family VCF files
 */

// Include modules
include { GLNEXUS_FAMILY } from '../modules/glnexus_family'

workflow FAMILY_CALLING {

    take:
    family_gvcfs    // channel: [fid, [barcodes], [gvcfs], [tbis]]

    main:

    // Run GLnexus for each family
    GLNEXUS_FAMILY(family_gvcfs)

    emit:
    family_vcfs = GLNEXUS_FAMILY.out.family_vcf
}