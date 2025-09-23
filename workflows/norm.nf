/*
 * VCF Normalization Workflow
 * Normalizes family VCF files using bcftools norm
 */

// Include modules
include { NORMALIZE } from '../modules/normalize_vcf'

workflow NORM {

    take:
    family_vcfs    // channel: [fid, vcf, tbi]

    main:

    // Run normalization for each family VCF
    NORMALIZE(family_vcfs)

    emit:
    normalized_vcfs = NORMALIZE.out.normalized_vcf
}