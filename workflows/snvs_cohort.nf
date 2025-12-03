/*
 * SNVs Cohort Workflow
 * Merges all family common_gt.bcf files into a cohort-level BCF
 * Concatenates all family DNM TSV files into a cohort-level TSV
 */

// Include modules
include { SNVS_COHORT_MERGE } from '../modules/snvs_cohort/snvs_cohort_merge'
include { MERGE_DNM } from '../modules/snvs_cohort/merge_dnm'

workflow SNVS_COHORT {

    take:
    filtered_common_bcfs    // channel: [fid, bcf, csi]
    dnm_files              // channel: [fid, bcf, csi, tsv]
    need_bcf_merge         // boolean: true if BCF merge is needed
    need_dnm_merge         // boolean: true if DNM merge is needed

    main:

    // Conditionally run BCF merge
    if (need_bcf_merge) {
        // Collect all BCF files and their indices for cohort merge
        bcf_files = filtered_common_bcfs
            .map { fid, bcf, csi -> bcf }
            .collect()
        
        csi_files = filtered_common_bcfs
            .map { fid, bcf, csi -> csi }
            .collect()

        // Run cohort merge
        SNVS_COHORT_MERGE(params.cohort_name, bcf_files, csi_files)
        cohort_bcf_output = SNVS_COHORT_MERGE.out.cohort_bcf
    } else {
        cohort_bcf_output = Channel.empty()
    }

    // Conditionally run DNM merge
    if (need_dnm_merge) {
        // Collect all DNM TSV files for concatenation
        dnm_tsv_files = dnm_files
            .map { fid, bcf, csi, tsv -> tsv }
            .collect()

        // Run DNM merge
        MERGE_DNM(params.cohort_name, dnm_tsv_files, params.vep_config_name)
        cohort_dnm_output = MERGE_DNM.out.cohort_dnm_tsv
    } else {
        cohort_dnm_output = Channel.empty()
    }

    emit:
    cohort_bcf = cohort_bcf_output
    cohort_dnm_tsv = cohort_dnm_output
}