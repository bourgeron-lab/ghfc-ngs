/*
 * Wombat Workflow
 * Integrated workflow for variant annotation and filtering:
 * 1. Annotate with gnomAD frequencies
 * 2. Filter into rare and common variants
 * 3. Filter common variants (keep only GT)
 * 4. Annotate rare variants with VEP
 */

// Include modules
include { W_GNOMAD_FREQ_ANNOT } from '../modules/w_gnomad_freq_annot'
include { W_GNOMAD_FREQ_FILTER } from '../modules/w_gnomad_freq_filter'
include { W_COMMON_FILTERS } from '../modules/w_common_filters'
include { W_VEP_ANNOTATION } from '../modules/w_vep_annotation'

workflow WOMBAT {

    take:
    normalized_bcfs    // channel: [fid, bcf, csi]

    main:

    // Step 1: Annotate with gnomAD frequencies
    gnomad_input = normalized_bcfs
        .map { fid, bcf, csi ->
            [fid, bcf, csi, params.gnomad_file]
        }

    W_GNOMAD_FREQ_ANNOT(gnomad_input)

    // Step 2: Filter into rare and common variants
    filter_input = W_GNOMAD_FREQ_ANNOT.out.annotated_bcfs
        .map { fid, bcf, csi ->
            [fid, bcf, csi, params.gnomad_filter_field, params.gnomad_filter_threshold]
        }

    W_GNOMAD_FREQ_FILTER(filter_input)

    // Step 3: Filter common variants (keep only GT field)
    W_COMMON_FILTERS(W_GNOMAD_FREQ_FILTER.out.common_bcfs)

    // Step 4: Annotate rare variants with VEP
    vep_input = W_GNOMAD_FREQ_FILTER.out.rare_bcfs
        .map { fid, bcf, csi ->
            [fid, bcf, csi, params.vep_config]
        }

    W_VEP_ANNOTATION(vep_input)

    emit:
    gnomad_bcfs = W_GNOMAD_FREQ_ANNOT.out.annotated_bcfs
    rare_bcfs = W_GNOMAD_FREQ_FILTER.out.rare_bcfs
    common_bcfs = W_GNOMAD_FREQ_FILTER.out.common_bcfs
    filtered_common_bcfs = W_COMMON_FILTERS.out.filtered_common_bcfs
    vep_annotated_bcfs = W_VEP_ANNOTATION.out.annotated_bcfs
}
