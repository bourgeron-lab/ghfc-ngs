/*
 * Wombat Workflow
 * Integrated workflow for variant annotation and filtering:
 * 1. Annotate with gnomAD frequencies
 * 2. Filter into rare and common variants
 * 3. Filter common variants (keep only GT)
 * 4. Annotate rare variants with VEP
 */

// Include modules
include { GNOMAD_FREQ_ANNOT } from '../modules/wombat/gnomad_freq_annot'
include { GNOMAD_FREQ_FILTER } from '../modules/wombat/gnomad_freq_filter'
include { COMMON_FILTERS } from '../modules/wombat/common_filters'
include { VEP_ANNOTATION } from '../modules/wombat/vep_annotation'
include { BCFTOOLS_ANNOTATE } from '../modules/wombat/bcftools_annotate'

workflow WOMBAT {

    take:
    normalized_bcfs    // channel: [fid, bcf, csi]

    main:

    // Step 1: Annotate with gnomAD frequencies
    gnomad_input = normalized_bcfs
        .map { fid, bcf, csi ->
            [fid, bcf, csi, params.gnomad_file]
        }

    GNOMAD_FREQ_ANNOT(gnomad_input)

    // Step 2: Filter into rare and common variants
    filter_input = GNOMAD_FREQ_ANNOT.out.annotated_bcfs
        .map { fid, bcf, csi ->
            [fid, bcf, csi, params.gnomad_filter_field, params.gnomad_filter_threshold]
        }

    GNOMAD_FREQ_FILTER(filter_input)

    // Step 3: Filter common variants (keep only GT field)
    COMMON_FILTERS(GNOMAD_FREQ_FILTER.out.common_bcfs)

    // Step 4: Annotate rare variants with VEP
    vep_input = GNOMAD_FREQ_FILTER.out.rare_vcfs
        .map { fid, vcf, tbi ->
            [fid, vcf, tbi, params.vep_config]
        }

    VEP_ANNOTATION(vep_input)

    // Step 5: Add additional annotations from BCF files
    bcftools_annotate_input = VEP_ANNOTATION.out.annotated_vcfs
        .map { fid, vcf, tbi ->
            [fid, vcf, tbi, params.wombat_annotation_path, params.wombat_annotation_list]
        }

    BCFTOOLS_ANNOTATE(bcftools_annotate_input)

    emit:
    gnomad_bcfs = GNOMAD_FREQ_ANNOT.out.annotated_bcfs
    rare_vcfs = GNOMAD_FREQ_FILTER.out.rare_vcfs
    common_bcfs = GNOMAD_FREQ_FILTER.out.common_bcfs
    filtered_common_bcfs = COMMON_FILTERS.out.filtered_common_bcfs
    vep_annotated_vcfs = VEP_ANNOTATION.out.annotated_vcfs
    fully_annotated_bcfs = BCFTOOLS_ANNOTATE.out.annotated_bcfs
}
