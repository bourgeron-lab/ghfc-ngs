/*
 * SNVs Common Filters Workflow
 * Filters common variants and keeps only GT format field
 */

// Include modules
include { snvs_common_filters } from '../modules/snvs_common_filters'

workflow SNVS_COMMON_FILTERS {

    take:
    common_vcfs    // channel: [fid, vcf, tbi]

    main:

    // Run common variants filtering
    snvs_common_filters(common_vcfs)

    emit:
    filtered_common_bcfs = snvs_common_filters.out.filtered_common_bcfs
}