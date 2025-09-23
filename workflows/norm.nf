/*
 * VCF Normalization Workflow
 * Normalizes family VCF files using bcftools norm
 */

// Include modules
include { NORM } from '../modules/normalize_vcf'

workflow NORM {

    take:
    family_vcfs    // channel: [fid, vcf, tbi]

    main:

    // Run normalization for each family VCF
    NORM(family_vcfs)

    emit:
    normalized_vcfs = NORM.out.normalized_vcf
}