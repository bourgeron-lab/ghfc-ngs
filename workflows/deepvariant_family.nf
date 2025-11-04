#!/usr/bin/env nextflow

/*
 * DeepVariant Family Workflow
 * Integrates GLnexus family calling, VCF normalization, and pedigree extraction
 */

// Include modules
include { GL_GLNEXUS_FAMILY } from '../modules/gl_glnexus_family'
include { GL_NORMALIZE } from '../modules/gl_normalize_vcf'
include { GL_FAMILIAL_PEDIGREE } from '../modules/gl_familial_pedigree'

workflow DEEPVARIANT_FAMILY {

    take:
    family_gvcfs    // channel: [fid, [barcodes], [gvcfs], [tbis]]
    pedigree_file   // path: cohort pedigree file

    main:

    // Run GLnexus for each family
    GL_GLNEXUS_FAMILY(family_gvcfs)

    // Normalize family VCF files
    GL_NORMALIZE(GL_GLNEXUS_FAMILY.out.family_vcf)

    // Extract family-specific pedigree files
    family_pedigree_input = GL_GLNEXUS_FAMILY.out.family_vcf
        .map { fid, vcf, tbi ->
            [fid, pedigree_file, params.cohort_name]
        }
    
    GL_FAMILIAL_PEDIGREE(family_pedigree_input)

    emit:
    family_vcfs = GL_GLNEXUS_FAMILY.out.family_vcf
    normalized_bcfs = GL_NORMALIZE.out.normalized_bcf
    family_pedigrees = GL_FAMILIAL_PEDIGREE.out.family_pedigree
}
