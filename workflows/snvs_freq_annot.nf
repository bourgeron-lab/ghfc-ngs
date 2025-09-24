/*
 * SNVs Frequency Annotation Workflow
 * Annotates normalized VCF files with gnomAD frequency data
 */

// Include modules
include { snvs_gnomad_freq_annot } from '../modules/snvs_gnomad_freq_annot'

workflow SNVS_FREQ_ANNOT {

    take:
    normalized_vcfs    // channel: [fid, vcf, tbi]

    main:

    // Prepare gnomAD annotation input channel - add gnomad_file to each tuple
    gnomad_input = normalized_vcfs
        .map { fid, vcf, tbi ->
            [fid, vcf, tbi, params.gnomad_file]
        }

    // Run gnomAD frequency annotation
    snvs_gnomad_freq_annot(gnomad_input)

    emit:
    annotated_vcfs = snvs_gnomad_freq_annot.out.annotated_vcfs
}