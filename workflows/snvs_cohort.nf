/*
 * SNVs Cohort Workflow
 * Merges all family common_gt.bcf files into a cohort-level BCF
 * Concatenates all family DNM TSV files into a cohort-level TSV
 */

// Include modules
include { SNVS_COHORT_MERGE } from '../modules/snvs_cohort/snvs_cohort_merge'
include { MERGE_WOMBAT } from '../modules/snvs_cohort/merge_wombat'

workflow SNVS_COHORT {

    take:
    filtered_common_bcfs    // channel: [fid, bcf, csi]
    wombat_files           // channel: [fid, wombat_config_name, tsv]
    need_bcf_merge         // boolean: true if BCF merge is needed
    need_wombat_merges     // map: [wombat_config_name: boolean] - whether merge is needed for each config

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
    
    // Conditionally run WOMBAT merges for each config
    if (need_wombat_merges && !need_wombat_merges.isEmpty()) {
        // Group wombat files by config name
        wombat_grouped = wombat_files
            .groupTuple(by: 1)  // Group by wombat_config_name
            .filter { fid_list, wombat_config_name, file_list ->
                // Only process configs that need merging
                need_wombat_merges[wombat_config_name] == true
            }
            .map { fid_list, wombat_config_name, file_list ->
                tuple(params.cohort_name, file_list, params.vep_config_name, wombat_config_name, "results")
            }
        
        MERGE_WOMBAT(
            wombat_grouped.map { it[0] },  // cohort_name
            wombat_grouped.map { it[1] },  // files
            wombat_grouped.map { it[2] },  // vep_config_name
            wombat_grouped.map { it[3] },  // wombat_config_name
            wombat_grouped.map { it[4] }   // output_name
        )
        
        cohort_wombat_output = MERGE_WOMBAT.out.cohort_wombat_tsv
    } else {
        cohort_wombat_output = Channel.empty()
    }

    emit:
    cohort_bcf = cohort_bcf_output
    cohort_wombat_tsv = cohort_wombat_output
}