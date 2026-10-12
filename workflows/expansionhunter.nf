/*
 * ExpansionHunter Workflow
 * Genotypes the variant catalog's short tandem repeats in each sample's CRAM
 * Merges each family's calls into one VCF and a genotype table, classed against published thresholds
 * Merges the family files into cohort-level files, with each allele's place among the cohort's founders
 *
 * Skips processing if output files already exist - only runs necessary modules
 */

// Include modules
include { EXPANSIONHUNTER_CALL } from '../modules/expansionhunter/expansionhunter_call'
include { EXPANSIONHUNTER_FAMILY } from '../modules/expansionhunter/expansionhunter_family'
include { EXPANSIONHUNTER_COHORT } from '../modules/expansionhunter/expansionhunter_cohort'

workflow EXPANSIONHUNTER {

    take:
    cram_files         // channel: [barcode, cram, crai] - existing and newly aligned CRAMs
    sex                // map: [barcode: 'male'|'female'|''] - from the pedigree
    family_members     // map: [barcode: fid]
    pedigree_file      // string: path to the cohort pedigree
    need_call          // collection: samples that need an ExpansionHunter call
    existing_calls     // collection: samples whose call is already on disk
    need_family        // collection: families whose VCF and table must be (re)built
    existing_families  // collection: families whose VCF and table are on disk and current
    need_cohort_merge  // boolean: whether the cohort files must be (re)built

    main:

    def catalog_name = params.expansionhunter_catalog_name
    def ref = file(params.ref)
    def ref_fai = file("${params.ref}.fai")
    def catalog = file(params.expansionhunter_catalog)
    // An empty list stages nothing, and the table is then written without classes
    def thresholds = params.expansionhunter_thresholds ? file(params.expansionhunter_thresholds) : []

    // Staged as process inputs rather than referenced through projectDir: the pipeline
    // directory is not mounted inside the task container.
    def family_table = file("${projectDir}/modules/expansionhunter/scripts/eh_family_table", checkIfExists: true)
    def cohort_merge = file("${projectDir}/modules/expansionhunter/scripts/eh_cohort_merge", checkIfExists: true)

    // A family's table holds every member that has, or gets in this run, a call; a member
    // with neither a call nor a CRAM is left out of it, as the plan already warned
    def expected_members = [:]
    family_members.each { barcode, fid ->
        if (fid in need_family && (barcode in need_call || barcode in existing_calls)) {
            expected_members.get(fid, []) << barcode
        }
    }

    // The cohort files hold every family with at least one member that has a call
    def expected_families = [] as Set
    expected_families.addAll(need_family)
    expected_families.addAll(existing_families)

    // =====================================================================
    // Step 1: Genotype each sample
    // =====================================================================
    EXPANSIONHUNTER_CALL(
        cram_files
            .filter { barcode, _cram, _crai -> barcode in need_call }
            .unique { barcode, _cram, _crai -> barcode }
            .map { barcode, cram, crai ->
                tuple(barcode, cram, crai, sex[barcode] ?: '', ref, ref_fai, catalog, catalog_name)
            }
    )

    // =====================================================================
    // Step 2: Merge and tabulate each family
    // =====================================================================
    if (need_family) {
        // Members called in an earlier run, looked up at their sharded path, never globbed
        existing_member_vcfs = channel.fromList(
            expected_members.collectMany { fid, barcodes ->
                barcodes.findAll { barcode -> barcode in existing_calls }.collect { barcode ->
                    tuple(barcode, file("${Sharding.getSampleDir(params.data, barcode)}/str/expansionhunter/${barcode}.${catalog_name}.vcf"))
                }
            })

        // Keyed with the family's expected size, so a family starts as soon as its last member
        // is called rather than when every sample of the cohort is. remainder keeps a family
        // that never fills up, so that the check below can name it instead of it vanishing.
        family_inputs = EXPANSIONHUNTER_CALL.out.vcf
            .mix(existing_member_vcfs)
            .filter { barcode, _vcf -> family_members[barcode] in need_family }
            .map { barcode, vcf ->
                def fid = family_members[barcode]
                tuple(groupKey(fid, expected_members[fid].size()), barcode, vcf)
            }
            .groupTuple(by: 0, remainder: true)
            .filter { fid, barcodes, _vcfs ->
                // Never publish a family table built from part of the family
                def missing = expected_members[fid.toString()] - barcodes
                if (missing) {
                    log.warn "Skipping ExpansionHunter family merge of ${fid}: no call for ${missing.join(', ')}"
                    return false
                }
                return true
            }
            .map { fid, _barcodes, vcfs ->
                tuple(fid.toString(), vcfs, file(pedigree_file), ref_fai, thresholds, catalog_name, family_table)
            }

        EXPANSIONHUNTER_FAMILY(family_inputs)
        new_family_outputs = EXPANSIONHUNTER_FAMILY.out.table
            .join(EXPANSIONHUNTER_FAMILY.out.vcf)
        new_family_vcfs = EXPANSIONHUNTER_FAMILY.out.vcf
        new_family_tables = EXPANSIONHUNTER_FAMILY.out.table
    } else {
        new_family_outputs = channel.empty()
        new_family_vcfs = channel.empty()
        new_family_tables = channel.empty()
    }

    // =====================================================================
    // Step 3: Merge the families into cohort files
    // =====================================================================
    if (need_cohort_merge) {
        // Families not rebuilt in this run contribute the files they already have
        existing_family_outputs = channel.fromList(
            existing_families.collect { fid ->
                def prefix = "${Sharding.getFamilyDir(params.data, fid)}/str/expansionhunter/${fid}.${catalog_name}"
                tuple(fid, file("${prefix}.genotypes.tsv"), file("${prefix}.vcf.gz"), file("${prefix}.vcf.gz.tbi"))
            })

        // Concat, then a single collect, so the merge waits for this run's own families. Never
        // publish cohort files built from an incomplete set of families: one whose merge
        // failed or was skipped blocks the cohort rather than dropping out of it unnoticed.
        cohort_input = new_family_outputs
            .concat(existing_family_outputs)
            .unique { fid, _table, _vcf, _tbi -> fid }
            .toList()
            .filter { rows ->
                def missing = expected_families - (rows.collect { fid, _table, _vcf, _tbi -> fid } as Set)
                if (missing) {
                    log.warn "Skipping ExpansionHunter cohort merge: no family files for ${missing.sort().take(20).join(', ')}${missing.size() > 20 ? ' and others' : ''}"
                    return false
                }
                return !rows.isEmpty()
            }
            .map { rows ->
                tuple(params.cohort_name,
                      rows.collect { _fid, table, _vcf, _tbi -> table },
                      rows.collect { _fid, _table, vcf, _tbi -> vcf },
                      rows.collect { _fid, _table, _vcf, tbi -> tbi },
                      catalog_name, cohort_merge)
            }

        EXPANSIONHUNTER_COHORT(cohort_input)
        cohort_vcf_output = EXPANSIONHUNTER_COHORT.out.vcf
        cohort_tables_output = EXPANSIONHUNTER_COHORT.out.tables
    } else {
        cohort_vcf_output = channel.empty()
        cohort_tables_output = channel.empty()
    }

    emit:
    sample_vcfs = EXPANSIONHUNTER_CALL.out.vcf
    family_vcfs = new_family_vcfs
    family_tables = new_family_tables
    cohort_vcf = cohort_vcf_output
    cohort_tables = cohort_tables_output
}
