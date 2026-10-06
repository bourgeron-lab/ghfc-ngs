/*
 * SNVs Cohort Workflow
 * Merges all family common_gt.bcf files into a cohort-level BCF
 * Concatenates all family DNM TSV files into a cohort-level TSV
 */

// Include modules
include { SNVS_COHORT_MERGE } from '../modules/snvs_cohort/snvs_cohort_merge'

workflow SNVS_COHORT {

    take:
    filtered_common_bcfs    // channel: [fid, bcf, csi]
    need_bcf_merge         // boolean: true if BCF merge is needed
    families               // list: every family in the pedigree
    annotation_needed      // list: families whose common filtered BCF this run should produce

    main:

    // Conditionally run BCF merge. The cohort wombat tables are the wombat_cohort step's.
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
    
    emit:
    cohort_bcf = cohort_bcf_output
}