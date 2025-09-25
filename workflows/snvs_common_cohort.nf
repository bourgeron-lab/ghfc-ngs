/*
 * SNVs Common Cohort Workflow
 * Merges all family common_gt.bcf files into a cohort-level BCF
 */

// Include modules
include { snvs_common_cohort } from '../modules/snvs_common_cohort'

workflow SNVS_COMMON_COHORT {

    take:
    filtered_common_bcfs    // channel: [fid, bcf, csi]

    main:

    // Collect all BCF files and their indices
    bcf_files = filtered_common_bcfs
        .map { fid, bcf, csi -> bcf }
        .collect()
    
    csi_files = filtered_common_bcfs
        .map { fid, bcf, csi -> csi }
        .collect()

    // Run cohort merge
    snvs_common_cohort(params.cohort_name, bcf_files, csi_files)

    emit:
    cohort_bcf = snvs_common_cohort.out.cohort_bcf
}