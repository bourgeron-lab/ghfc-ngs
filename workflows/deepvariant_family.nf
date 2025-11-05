#!/usr/bin/env nextflow

/*
 * DeepVariant Family Workflow
 * Integrates GLnexus family calling, VCF normalization, and pedigree extraction
 */

// Include modules
include { GLNEXUS_FAMILY } from '../modules/deepvariant_family/glnexus_family'
include { NORMALIZE } from '../modules/deepvariant_family/normalize_vcf'
include { FAMILIAL_PEDIGREE } from '../modules/deepvariant_family/familial_pedigree'

workflow DEEPVARIANT_FAMILY {

    take:
    family_gvcfs    // channel: [fid, [barcodes], [gvcfs], [tbis]]
    pedigree_file   // path: cohort pedigree file

    main:

    // Run GLnexus for each family
    GLNEXUS_FAMILY(family_gvcfs)

    // Normalize family VCF files
    NORMALIZE(GLNEXUS_FAMILY.out.family_vcf)

    // Extract family-specific pedigree files
    family_pedigree_input = GLNEXUS_FAMILY.out.family_vcf
        .map { fid, vcf, tbi ->
            [fid, pedigree_file, params.cohort_name]
        }
    
    FAMILIAL_PEDIGREE(family_pedigree_input)

    emit:
    family_vcfs = GLNEXUS_FAMILY.out.family_vcf
    normalized_bcfs = NORMALIZE.out.normalized_bcf
    family_pedigrees = FAMILIAL_PEDIGREE.out.family_pedigree
}
