/*
 * SNVs Frequency Filtering Workflow
 * Filters gnomAD-annotated VCF files based on frequency thresholds
 */

// Include modules
include { snvs_gnomad_freq_filter } from '../modules/snvs_gnomad_freq_filter'

workflow SNVS_FREQ_FILTER {

    take:
    gnomad_vcfs    // channel: [fid, vcf, tbi]

    main:

    // Prepare filtering input channel - add filter parameters to each tuple
    filter_input = gnomad_vcfs
        .map { fid, vcf, tbi ->
            [fid, vcf, tbi, params.gnomad_filter_field, params.gnomad_filter_threshold]
        }

    // Run gnomAD frequency filtering
    snvs_gnomad_freq_filter(filter_input)

    emit:
    filtered_vcfs = snvs_gnomad_freq_filter.out.filtered_vcfs
}