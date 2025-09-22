/*
 * VEP Annotation Workflow
 * Annotates family VCF files using Ensembl VEP
 */

// Include modules
include { VEP_ANNOTATION } from '../modules/vep_annotation'

workflow VEP_ANNOTATION_WORKFLOW {

    take:
    family_vcfs    // channel: [fid, vcf, tbi]

    main:

    // Run VEP annotation for each family VCF
    VEP_ANNOTATION(family_vcfs)

    emit:
    annotated_vcfs = VEP_ANNOTATION.out.annotated_vcf
}