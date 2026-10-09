/*
 * ROH / Inbreeding Workflow
 * Calls runs of homozygosity on each family's ancestry panel genotypes
 * Scores every member's inbreeding against the ancestry-pgs reference
 * Concatenates the family tables into cohort-level tables
 *
 * Skips processing if output files already exist - only runs necessary modules
 */

// Include modules
include { ROH_CALL } from '../modules/roh/roh_call'
include { ROH_COHORT_MERGE } from '../modules/roh/roh_cohort_merge'

workflow ROH {

    take:
    family_inputs      // channel: [fid, panel_bcf, panel_csi, ancestry_tsv] - one per family needing a call
    pedigree_file      // string: path to the cohort pedigree
    families           // collection: family IDs in the current pedigree - scopes the cohort merge
    need_call          // collection: families that need ROH calling
    need_cohort_merge  // boolean: whether the cohort table merge is needed
    gene_beds          // list: [symbol BED, ID BED], both empty when no gene annotation is configured

    main:

    def panel_name = params.ancestry_panel_name
    def reference = params.ancestry_reference
    def (genes_bed, gene_ids_bed) = gene_beds

    // Staged as a process input rather than referenced through projectDir: the
    // pipeline directory is not mounted inside the task container.
    def roh_cohort_merge = file("${projectDir}/modules/roh/scripts/roh_cohort_merge", checkIfExists: true)

    // =====================================================================
    // Step 1: Call and score each family
    // =====================================================================
    ROH_CALL(
        family_inputs
            .filter { fid, _bcf, _csi, _anc -> fid in need_call }
            .unique { fid, _bcf, _csi, _anc -> fid }
            .map { fid, bcf, csi, anc ->
                tuple(fid, bcf, csi, anc, file(pedigree_file), reference, panel_name,
                      genes_bed, gene_ids_bed, params.roh_min_gene_roh_mb)
            }
    )

    // =====================================================================
    // Step 2: Concatenate the family tables into cohort tables
    // =====================================================================
    if (need_cohort_merge) {
        def kinds = ['froh', 'roh', 'roh_genes', 'kinship']

        new_family_tables = ROH_CALL.out.froh
            .mix(ROH_CALL.out.roh, ROH_CALL.out.roh_genes, ROH_CALL.out.kinship)

        // Families not called in this run contribute the tables they already have
        existing_family_tables = channel.fromList(
            families
                .findAll { fid -> !(fid in need_call) }
                .collectMany { fid ->
                    kinds.collect { kind ->
                        file("${Sharding.getFamilyDir(params.data, fid)}/roh/${fid}.${panel_name}.${kind}.tsv")
                    }.findAll { tsv -> tsv.exists() }
                     .collect { tsv -> tuple(fid, tsv) }
                })

        // Concat, then a single collect, so the merge waits for this run's own calls.
        // Never publish a cohort table built from an incomplete set of families: a
        // family that should have a per-sample table and has none blocks the merge
        // rather than dropping out of the cohort unnoticed.
        cohort_merge_input = new_family_tables
            .concat(existing_family_tables)
            .unique { _fid, tsv -> tsv.name }
            .toList()
            .filter { rows ->
                def have = rows.findAll { _fid, tsv -> tsv.name.endsWith('.froh.tsv') }
                               .collect { fid, _tsv -> fid } as Set
                def expected = families.findAll { fid ->
                    fid in need_call ||
                    file("${Sharding.getFamilyDir(params.data, fid)}/roh/${fid}.${panel_name}.froh.tsv").exists()
                }
                def missing = expected - have
                if (missing) {
                    log.warn "Skipping cohort ROH merge: ROH tables missing for ${missing.join(', ')}"
                    return false
                }
                return !have.isEmpty()
            }
            .map { rows ->
                tuple(params.cohort_name, rows.collect { _fid, tsv -> tsv }, panel_name, roh_cohort_merge)
            }

        ROH_COHORT_MERGE(cohort_merge_input)
        cohort_tables_output = ROH_COHORT_MERGE.out.cohort_tables
    } else {
        cohort_tables_output = channel.empty()
    }

    emit:
    family_froh = ROH_CALL.out.froh
    family_roh = ROH_CALL.out.roh
    family_roh_genes = ROH_CALL.out.roh_genes
    cohort_tables = cohort_tables_output
}
