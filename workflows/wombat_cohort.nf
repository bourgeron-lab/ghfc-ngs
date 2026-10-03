/*
 * Wombat Cohort Workflow
 * Concatenates every family's wombat table into one cohort table per wombat config
 *
 * A step of its own, apart from snvs_cohort, so that a cohort can have its wombat tables
 * without the common-variant BCF merge - for a cohort of tens of thousands of exome families,
 * that merge is the expensive part, and the wombat tables are what gets read.
 */

// Include modules
include { MERGE_WOMBAT } from '../modules/wombat_cohort/merge_wombat'

workflow WOMBAT_COHORT {

    take:
    wombat_files           // channel: [fid, wombat_config_name, tsv]
    need_merges            // map: [wombat_config_name: boolean] - whether merge is needed for each config
    families               // list: every family in the pedigree
    wombat_needed          // list: families whose wombat tables this run should produce

    main:

    // Conditionally run the merge of each config
    if (need_merges && !need_merges.isEmpty()) {
        // Group wombat files by config name
        wombat_grouped = wombat_files
            .groupTuple(by: 1)  // Group by wombat_config_name
            .filter { fid_list, wombat_config_name, file_list ->
                // Only process configs that need merging
                need_merges[wombat_config_name] == true
            }
            .filter { fid_list, wombat_config_name, _file_list ->
                // Never publish a cohort table built from an incomplete set of families.
                // A family that produced nothing in this run must block the merge rather
                // than drop out of the cohort table unnoticed.
                def expected = families.findAll { fid ->
                    fid in wombat_needed ||
                    file("${Sharding.getFamilyDir(params.data, fid)}/wombat/${fid}.rare.${params.vep_config_name}.annotated.${wombat_config_name}.tsv").exists()
                }
                def missing = expected - fid_list
                if (missing) {
                    log.warn "Skipping cohort ${wombat_config_name} merge: wombat table missing for ${missing.join(', ')}"
                    return false
                }
                return true
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
    cohort_wombat_tsv = cohort_wombat_output
}
