/*
 * Ancestry and PGS Workflow
 * Extracts the reference panel sites from each sample's DeepVariant gVCF
 * Merges the per-sample panel genotypes per family
 * Projects each family onto the reference panel and scores the PGS catalog
 * Concatenates the family tables into cohort-level tables
 *
 * Skips processing if output files already exist - only runs necessary modules
 */

// Include modules
include { PANEL_SITES } from '../modules/ancestry/panel_sites'
include { PANEL_EXTRACT } from '../modules/ancestry/panel_extract'
include { PANEL_MERGE_FAMILY } from '../modules/ancestry/panel_merge_family'
include { ANCESTRY_PCS } from '../modules/ancestry/ancestry_pcs'
include { ANCESTRY_ADMIXTURE } from '../modules/ancestry/ancestry_admixture'
include { ANCESTRY_PGS_RAW } from '../modules/ancestry/ancestry_pgs_raw'
include { ANCESTRY_PGS_APPLY } from '../modules/ancestry/ancestry_pgs_apply'
include { ANCESTRY_COHORT_MERGE } from '../modules/ancestry/merge_cohort_tables'

workflow ANCESTRY {

    take:
    gvcfs              // channel: [barcode, gvcf, tbi]
    family_members     // map: [barcode: fid] - mapping of barcodes to family IDs
    families           // collection: family IDs in the current pedigree - scopes the cohort merge
    need_extract       // list of barcodes that need panel extraction
    need_family_merge  // map: [fid: boolean] - whether the family panel merge is needed
    need_family_score  // map: [fid: boolean] - whether family scoring is needed
    need_cohort_merge  // boolean: whether the cohort table merge is needed

    main:

    // The per-family tables that get concatenated into cohort tables. These are the
    // table-name suffixes ancestry-pgs writes; none is a suffix of another, which the
    // merge module's glob relies on.
    def table_kinds = ['pcs', 'ancestry', 'Q', 'pgs_raw', 'pgs_adjusted', 'pgs_zscore']

    def panel_name = params.ancestry_panel_name
    def reference = params.ancestry_reference
    def catalog = params.ancestry_catalog

    // Expected sample count per family, so a partial family is never merged
    def family_sizes = family_members.values().countBy { fid -> fid }

    // =====================================================================
    // Step 1: Verify the bundle and its fitted model, build the union site list
    //         (once per run)
    // =====================================================================
    PANEL_SITES(channel.of(tuple(params.cohort_name, reference, catalog, panel_name)))
    panel_sites = PANEL_SITES.out.panel_sites

    // =====================================================================
    // Step 2: Extract the panel sites from each sample's gVCF
    // =====================================================================
    extract_input = gvcfs
        .filter { barcode, _gvcf, _tbi -> barcode in need_extract }
        .combine(panel_sites)
        .map { barcode, gvcf, tbi, _pname, sites_file, regions_file ->
            tuple(barcode, gvcf, tbi, sites_file, regions_file, reference, panel_name)
        }

    PANEL_EXTRACT(extract_input)

    // Existing per-sample panel genotypes already on disk, restricted to the samples
    // in the current pedigree. Without that filter the glob picks up every sample
    // ever processed under params.data, including other cohorts'.
    existing_sample_panels = channel
        .fromPath("${params.data}/samples/*/*/*/ancestry/*.panel_gt.${panel_name}.bcf")
        .map { bcf ->
            def barcode = bcf.name - ".panel_gt.${panel_name}.bcf"
            tuple(barcode, bcf, file("${bcf}.csi"))
        }
        .filter { barcode, _bcf, csi ->
            family_members.containsKey(barcode) && csi.exists()
        }

    // A sample extracted in this run also matches the disk glob once published, so
    // dedupe by barcode - staging both copies would make the family merge count the
    // same sample twice and fail its sample-count check.
    all_sample_panels = PANEL_EXTRACT.out.panel_gt
        .mix(existing_sample_panels)
        .unique { barcode, _bcf, _csi -> barcode }

    // =====================================================================
    // Step 3: Merge the per-sample panel genotypes per family
    // =====================================================================
    family_merge_input = all_sample_panels
        .map { barcode, bcf, csi -> tuple(family_members[barcode], barcode, bcf, csi) }
        .groupTuple(by: 0)
        .filter { fid, barcodes, _bcfs, _csis ->
            if (need_family_merge[fid] != true) {
                return false
            }
            // Never publish a family panel built from part of the family: a missing
            // member would silently drop out of every downstream table.
            def expected = family_sizes[fid] ?: 0
            if (barcodes.size() != expected) {
                log.warn "Skipping ancestry panel merge for family ${fid}: ${barcodes.size()} of ${expected} samples have panel genotypes"
                return false
            }
            return true
        }
        .map { fid, _barcodes, bcfs, csis ->
            tuple(fid, bcfs, csis, panel_name, family_sizes[fid])
        }

    PANEL_MERGE_FAMILY(family_merge_input)

    existing_family_panels = channel
        .fromPath("${params.data}/families/*/*/*/ancestry/*.panel_gt.${panel_name}.bcf")
        .map { bcf ->
            def fid = bcf.name - ".panel_gt.${panel_name}.bcf"
            tuple(fid, bcf, file("${bcf}.csi"))
        }
        .filter { fid, _bcf, csi -> fid in families && csi.exists() }

    all_family_panels = PANEL_MERGE_FAMILY.out.family_panel_gt
        .mix(existing_family_panels)
        .unique { fid, _bcf, _csi -> fid }

    // =====================================================================
    // Step 4: Score each family against the reference panel
    // =====================================================================
    families_to_score = all_family_panels
        .filter { fid, _bcf, _csi -> need_family_score[fid] == true }

    ANCESTRY_PCS(
        families_to_score.map { fid, bcf, csi ->
            tuple(fid, bcf, csi, reference, panel_name)
        }
    )

    ANCESTRY_ADMIXTURE(
        families_to_score.map { fid, bcf, csi ->
            tuple(fid, bcf, csi, reference, panel_name)
        }
    )

    ANCESTRY_PGS_RAW(
        families_to_score.map { fid, bcf, csi ->
            tuple(fid, bcf, csi, reference, catalog, panel_name)
        }
    )

    // The adjustment must use the PCs from ANCESTRY_PCS and no other source:
    // projected components are on a different scale from in-sample ones, and the
    // fitted model is calibrated against the projected kind.
    ANCESTRY_PGS_APPLY(
        ANCESTRY_PGS_RAW.out.pgs_raw
            .join(ANCESTRY_PCS.out.pcs)
            .map { fid, pgs_raw, pcs ->
                tuple(fid, pgs_raw, pcs, reference, catalog, panel_name)
            }
    )

    // =====================================================================
    // Step 5: Concatenate the family tables into cohort tables
    // =====================================================================
    new_family_tables = ANCESTRY_PCS.out.pcs.map { fid, tsv -> tuple(fid, 'pcs', tsv) }
        .mix(ANCESTRY_PCS.out.ancestry.map { fid, tsv -> tuple(fid, 'ancestry', tsv) })
        .mix(ANCESTRY_ADMIXTURE.out.admixture.map { fid, tsv -> tuple(fid, 'Q', tsv) })
        .mix(ANCESTRY_PGS_RAW.out.pgs_raw.map { fid, tsv -> tuple(fid, 'pgs_raw', tsv) })
        .mix(ANCESTRY_PGS_APPLY.out.pgs_adjusted.map { fid, tsv -> tuple(fid, 'pgs_adjusted', tsv) })
        .mix(ANCESTRY_PGS_APPLY.out.pgs_zscore.map { fid, tsv -> tuple(fid, 'pgs_zscore', tsv) })

    if (need_cohort_merge) {
        existing_family_tables = channel.fromList(table_kinds)
            .flatMap { kind ->
                families.collect { fid ->
                    def tsv = file("${Sharding.getFamilyDir(params.data, fid)}/ancestry/${fid}.${panel_name}.${kind}.tsv")
                    tsv.exists() ? tuple(fid, kind, tsv) : null
                }.findAll { row -> row != null }
            }

        // Concat rather than collect-then-mix so the barrier genuinely waits for this
        // run's own scoring before the cohort tables are written. Dedupe in case a
        // freshly scored family has already been published to disk.
        cohort_merge_input = new_family_tables
            .concat(existing_family_tables)
            .unique { fid, kind, _tsv -> "${fid}:${kind}" }
            .groupTuple(by: 1)
            .filter { fids, kind, _tsvs ->
                // Never publish a cohort table built from an incomplete set of families.
                //
                // A family that could not be re-scored in this run still contributes the
                // table it already has, and that is correct: the panel label pins the
                // bundle and the thresholds, so an older table under the same label was
                // computed the same way. A family with no table at all is missing here
                // and blocks the merge rather than dropping out of the cohort unnoticed.
                def expected = families.findAll { fid ->
                    need_family_score[fid] == true ||
                    file("${Sharding.getFamilyDir(params.data, fid)}/ancestry/${fid}.${panel_name}.${kind}.tsv").exists()
                }
                def missing = expected - fids
                if (missing) {
                    log.warn "Skipping cohort ${kind} merge: ${kind} table missing for ${missing.join(', ')}"
                    return false
                }
                return true
            }
            .map { _fids, kind, tsvs ->
                tuple(params.cohort_name, kind, tsvs, panel_name)
            }

        ANCESTRY_COHORT_MERGE(cohort_merge_input)
        cohort_tables_output = ANCESTRY_COHORT_MERGE.out.cohort_table
    } else {
        cohort_tables_output = channel.empty()
    }

    emit:
    panel_sites = panel_sites
    sample_panel_gt = all_sample_panels
    sample_panel_stats = PANEL_EXTRACT.out.panel_stats
    family_panel_gt = all_family_panels
    family_pcs = ANCESTRY_PCS.out.pcs
    family_ancestry = ANCESTRY_PCS.out.ancestry
    family_admixture = ANCESTRY_ADMIXTURE.out.admixture
    family_pgs_zscore = ANCESTRY_PGS_APPLY.out.pgs_zscore
    cohort_tables = cohort_tables_output
}
