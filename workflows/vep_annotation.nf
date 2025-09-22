/*
 * Simple VEP Annotation Workflow
 * Annotates family VCF files using Ensembl VEP
 */

// Include the VEP annotation module
include { VEP_ANNOTATION } from '../modules/vep_annotation'

workflow VEP_ANNOTATION {

    take:
    family_vcfs    // channel: [fid, vcf, tbi]

    main:

    // Prepare VEP input channel - add vep_config to each tuple
    vep_input = family_vcfs
        .map { fid, vcf, tbi ->
            [fid, vcf, tbi, params.vep_config]
        }

    // Run VEP annotation
    VEP_ANNOTATION(vep_input)

    emit:
    annotated_vcfs = VEP_ANNOTATION.out.annotated_vcfs
}