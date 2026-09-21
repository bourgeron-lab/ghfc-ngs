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
    families               // list: every family in the pedigree
    annotation_needed      // list: families whose common filtered BCF this run should produce
    wombat_needed          // list: families whose wombat tables this run should produce

    main:

    // Conditionally run BCF merge
    if (need_bcf_merge) {
        // Collect the triples once so the barrier below filters the BCFs and their indices
        // consistently - collecting each of them independently cannot be filtered as a unit
        cohort_bcf_input = filtered_common_bcfs
            .toList()
            .filter { rows ->
                // toList emits [] on an empty channel where collect emitted nothing at all,
                // so keep the old behaviour of not running the merge with no inputs
                if (!rows) {
                    log.warn "Skipping cohort common variant merge: no family common filtered BCFs are available"
                    return false
                }
                // Never publish a cohort BCF built from an incomplete set of families
                def expected = families.findAll { fid ->
                    fid in annotation_needed ||
                    file("${Sharding.getFamilyDir(params.data, fid)}/vcfs/${fid}.common_gt.bcf").exists()
                }
                def missing = expected - rows.collect { fid, _bcf, _csi -> fid }
                if (missing) {
                    log.warn "Skipping cohort common variant merge: common filtered BCF missing for ${missing.join(', ')}"
                    return false
                }
                return true
            }

        // Run cohort merge
        SNVS_COHORT_MERGE(
            params.cohort_name,
            cohort_bcf_input.map { rows -> rows.collect { _fid, bcf, _csi -> bcf } },
            cohort_bcf_input.map { rows -> rows.collect { _fid, _bcf, csi -> csi } }
        )
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
    cohort_bcf = cohort_bcf_output
    cohort_wombat_tsv = cohort_wombat_output
}