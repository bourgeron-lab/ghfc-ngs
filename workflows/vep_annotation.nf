/*
 * Simple VEP Annotation Workflow
 * Annotates family VCF files using Ensembl VEP
 */

// Include the VEP annotation module
include { runVEP } from '../modules/vep_annotation'

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
    runVEP(vep_input)

    emit:
    annotated_vcfs = runVEP.out.annotated_vcfs
}