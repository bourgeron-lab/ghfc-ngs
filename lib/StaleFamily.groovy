/**
 * The rules behind --clean_stale_families: which files a family's joint call owns, and whether
 * a drifted family can actually be re-called from what is on disk and from the steps requested.
 *
 * Deliberately free of any Nextflow dependency - no `log`, no `params`, no `file()` - so that
 * every rule in here can be exercised with a plain `groovy -cp lib` script. That matters more
 * for this class than for any other in the repo: it is the only code in the pipeline that
 * decides what to destroy, and the caller in main.nf does the destroying.
 *
 * See the STALE FAMILY OUTPUTS handling in main.nf for how these lists are used.
 */
class StaleFamily {

    /**
     * Every published file a family owns that is derived from its joint call, in the order a
     * human wants to read them.
     *
     * Enumerated rather than globbed, on purpose: a glob over the family directory would also
     * catch files this pipeline never wrote, and the whole safety argument for deleting
     * anything is that we only remove what we know how to rebuild.
     *
     * Three groups are deliberately absent:
     *
     *  - `ancestry/` - it self-heals. A member with no gVCF has no panel BCF either, so it is
     *    already in the plan's ancestry `need_extract`, which forces `need_family_merge` and
     *    `need_family_score` for the family; and ANCESTRY refuses to publish a family panel
     *    built from part of a family. Deleting these would be pure loss.
     *  - `*.dnm.*` - DNM_EXTRACTION is included by no workflow, so nothing would regenerate
     *    them. A file left stale is recoverable; a file deleted that nothing rebuilds is not.
     *  - `extractor/` - the plan does no completeness check for it, so deleting would not
     *    cause a rebuild either.
     *
     * The annotation group is listed in full rather than just its final BCF because ANNOTATION
     * re-checks each of its stages independently: leaving `rare.vcf.gz` in place would have VEP
     * rebuilt from the stale call set, which is the exact bug this is meant to fix.
     */
    static List<String> familyOutputs(String data, String fid, String vepConfigName,
                                      List<String> wombatConfigNames, boolean includeWisecondorx) {
        def fam_dir = Sharding.getFamilyDir(data, fid)
        def vcfs = "${fam_dir}/vcfs"
        def paths = []

        // deepvariant_family - the three files whose presence reads as "this family is called"
        paths << "${vcfs}/${fid}.norm.bcf".toString()
        paths << "${vcfs}/${fid}.norm.bcf.csi".toString()
        paths << "${fam_dir}/${fid}.pedigree.tsv".toString()

        // annotation - all ten gate files, plus the gnomAD intermediate if it was left behind
        paths << "${vcfs}/${fid}.rare.vcf.gz".toString()
        paths << "${vcfs}/${fid}.rare.vcf.gz.tbi".toString()
        paths << "${vcfs}/${fid}.common.bcf".toString()
        paths << "${vcfs}/${fid}.common.bcf.csi".toString()
        paths << "${vcfs}/${fid}.common_gt.bcf".toString()
        paths << "${vcfs}/${fid}.common_gt.bcf.csi".toString()
        paths << "${vcfs}/${fid}.rare.${vepConfigName}.vcf.gz".toString()
        paths << "${vcfs}/${fid}.rare.${vepConfigName}.vcf.gz.tbi".toString()
        paths << "${vcfs}/${fid}.rare.${vepConfigName}.annotated.bcf".toString()
        paths << "${vcfs}/${fid}.rare.${vepConfigName}.annotated.bcf.csi".toString()
        paths << "${vcfs}/${fid}.gnomad.bcf".toString()
        paths << "${vcfs}/${fid}.gnomad.bcf.csi".toString()

        // wombat - the parquet, and one results TSV per configured wombat config
        paths << "${fam_dir}/wombat/${fid}.rare.${vepConfigName}.annotated.parquet".toString()
        configNames(wombatConfigNames).each { cfg ->
            paths << "${fam_dir}/wombat/${fid}.rare.${vepConfigName}.annotated.${cfg}.tsv".toString()
        }

        // wisecondorx - membership-dependent like the rest, but only when the step is in play
        if (includeWisecondorx) {
            paths << "${fam_dir}/svs/wisecondorx/${fid}_aberrations.bed".toString()
            paths << "${fam_dir}/svs/wisecondorx/${fid}_aberrations.annotated.bed".toString()
        }

        return paths
    }

    /**
     * The cohort-level merges that a family clean makes stale. The plan only checks that these
     * exist, so nothing would rebuild them on their own - they would keep data from a call set
     * that no longer exists, indefinitely.
     *
     * The cohort *ancestry* tables are not here: their `need_cohort_merge` is driven by the
     * per-family scores, so they already rebuild themselves.
     */
    static List<String> cohortOutputs(String data, String cohortName, String vepConfigName,
                                      List<String> wombatConfigNames, boolean includeWisecondorx) {
        def cohort_dir = "${data}/cohorts/${cohortName}"
        def paths = []

        paths << "${cohort_dir}/vcfs/${cohortName}.common_gt.bcf".toString()
        paths << "${cohort_dir}/vcfs/${cohortName}.common_gt.bcf.csi".toString()

        configNames(wombatConfigNames).each { cfg ->
            paths << "${cohort_dir}/wombat/${cohortName}.rare.${vepConfigName}.annotated.${cfg}.results.tsv".toString()
        }

        if (includeWisecondorx) {
            paths << "${cohort_dir}/svs/wisecondorx/${cohortName}_aberrations.bed".toString()
        }

        return paths
    }

    /**
     * Can this family's missing members actually be re-called?
     *
     * Deleting a family's outputs when the answer is no leaves the cohort strictly worse off
     * than the warning it replaces: the call set is gone, nothing schedules the missing member,
     * and the run aborts in validateStepsAvailability having already destroyed the evidence.
     *
     * Two independent halves. The data half asks whether each missing member can reach a gVCF:
     * one that already has a CRAM only needs variant calling, one without needs an alignment
     * input to resolve. The steps half asks whether this run was even asked to do that work -
     * a family the requested steps cannot rebuild is a guaranteed own-goal.
     *
     * Returns [recoverable: boolean, reasons: List<String>, needs_alignment: List<String>].
     */
    static Map recovery(Collection missing, Set haveCram, Set resolvable, Map indexMissing,
                        Collection steps) {
        def reasons = []
        def needs_alignment = []

        (missing ?: []).each { barcode ->
            if (haveCram?.contains(barcode)) {
                // A CRAM but no gVCF: variant calling alone gets this member back
                return
            }
            needs_alignment << barcode
            if (resolvable?.contains(barcode)) {
                return
            }
            // The data is there but unusable is a different problem with a different fix, so
            // it gets its own sentence rather than being folded into "no input data"
            def unindexed = indexMissing?.get(barcode)
            reasons << (unindexed
                ? "input CRAM present but unindexed for ${barcode}: ${unindexed.path} - run 'samtools index' on it".toString()
                : "no input data on disk for ${barcode}".toString())
        }

        def requested = (steps ?: []) as Set
        if (!requested.contains('deepvariant_family')) {
            reasons << "the deepvariant_family step was not requested, so the family could not be re-called".toString()
        }
        if (!requested.contains('deepvariant_sample')) {
            reasons << "the deepvariant_sample step was not requested, so the missing members could not be called".toString()
        }
        if (needs_alignment && !requested.contains('alignment')) {
            reasons << "the alignment step was not requested, but ${needs_alignment.join(', ')} would need aligning".toString()
        }

        return [recoverable: reasons.isEmpty(), reasons: reasons, needs_alignment: needs_alignment]
    }

    /** Wombat config file names minus their YAML suffix, the same way the plan derives them. */
    static List<String> configNames(List<String> wombatConfigNames) {
        if (!wombatConfigNames) return []
        return wombatConfigNames.collect { it?.toString()?.replaceAll(/\.ya?ml$/, '') }.findAll { it }
    }
}
