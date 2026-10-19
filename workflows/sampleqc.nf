/*
 * Sample QC Workflow
 * Extracts each sample's allele counts at somalier's sites from its CRAM
 * Checks each family's relatedness against its pedigree, and every member's sex and
 * karyotype, with the cnv_coverage JSONs when there are and the family BCF's Mendelian counts
 * Gathers the family tables for the cohort, with the relatedness between families
 *
 * Rebuilds a family when its inputs snapshot - its pedigree rows and the evidence used - is not
 * the one on disk, and the cohort when any family is rebuilt or its snapshot changes
 */

// Include modules
include { SOMALIER_EXTRACT } from '../modules/sampleqc/somalier_extract'
include { SAMPLEQC_FAMILY } from '../modules/sampleqc/sampleqc_family'
include { SAMPLEQC_COHORT } from '../modules/sampleqc/sampleqc_cohort'

workflow SAMPLEQC {

    take:
    cram_files        // channel: [barcode, cram, crai] - existing and newly aligned CRAMs
    new_coverage      // channel: [barcode, bw, json] - coverage the cnv_coverage step makes in this run
    norm_bcfs         // channel: [fid, bcf, csi] - existing and newly called family BCFs
    family_members    // map: [barcode: fid]
    qc                // map: the plan's sampleqc entry; see createAnalysisPlan

    main:

    def label = params.sampleqc_label
    def ref = file(params.ref)
    def ref_fai = file("${params.ref}.fai")
    def sites = file(params.somalier_sites)
    def cov_tag = "cov_b${params.cnv_cov_bin}q${params.cnv_cov_mapq}"

    // Staged as a process input rather than referenced through projectDir: the pipeline
    // directory is not mounted inside the task container.
    def classify = file("${projectDir}/modules/sampleqc/scripts/sampleqc_classify", checkIfExists: true)

    // Read from their sharded paths, never globbed
    def extract_file = { barcode ->
        file("${Sharding.getSampleDir(params.data, barcode)}/qc/somalier/${barcode}.${label}.somalier")
    }
    def coverage_file = { barcode ->
        file("${Sharding.getSampleDir(params.data, barcode)}/coverage/${barcode}.${cov_tag}.json")
    }
    def published_tables = { fid ->
        def prefix = "${Sharding.getFamilyDir(params.data, fid)}/qc/${fid}.${label}"
        ['samples', 'pairs', 'flags'].collect { kind -> file("${prefix}.${kind}.tsv") }
    }

    // =====================================================================
    // Step 1: Extract each sample
    // =====================================================================
    SOMALIER_EXTRACT(
        cram_files
            .filter { barcode, _cram, _crai -> barcode in qc.needed }
            .unique { barcode, _cram, _crai -> barcode }
            .map { barcode, cram, crai -> tuple(barcode, cram, crai, ref, ref_fai, sites, label) }
    )
    new_extracts = SOMALIER_EXTRACT.out.extract

    // =====================================================================
    // Step 2: Check each family
    // =====================================================================
    if (qc.need_family) {
        // Every input of a family, tagged with its kind: one extract and one coverage JSON per
        // member the snapshot says has them, and the family BCF when it says so
        def expected = [:]
        qc.need_family.each { fid ->
            expected[fid] = qc.extract_members[fid].size() + qc.coverage_members[fid].size() +
                (fid in qc.mendelian ? 1 : 0)
        }
        def from_disk = []
        qc.need_family.each { fid ->
            qc.extract_members[fid].findAll { barcode -> !(barcode in qc.needed) }.each { barcode ->
                from_disk << tuple(fid, 'extract', extract_file(barcode))
            }
            qc.coverage_members[fid].findAll { barcode -> !(barcode in qc.coverage_from_run) }.each { barcode ->
                from_disk << tuple(fid, 'coverage', coverage_file(barcode))
            }
        }

        family_items = channel.fromList(from_disk)
            .mix(new_extracts
                .filter { barcode, _extract -> family_members[barcode] in qc.need_family }
                .map { barcode, extract -> tuple(family_members[barcode], 'extract', extract) })
            .mix(new_coverage
                .filter { barcode, _bw, _json -> barcode in qc.coverage_from_run &&
                          family_members[barcode] in qc.need_family }
                .map { barcode, _bw, json -> tuple(family_members[barcode], 'coverage', json) })
            .mix(norm_bcfs
                .filter { fid, _bcf, _csi -> fid in qc.need_family && fid in qc.mendelian }
                .unique { fid, _bcf, _csi -> fid }
                .map { fid, bcf, csi -> tuple(fid, 'norm', [bcf, csi]) })

        // Keyed with the family's expected count, so a family starts as soon as its last input
        // is made. remainder keeps one that never fills up, so that the check below names it.
        family_inputs = family_items
            .map { fid, kind, payload -> tuple(groupKey(fid, expected[fid]), kind, payload) }
            .groupTuple(by: 0, remainder: true)
            .filter { fid, kinds, _payloads ->
                // Never check a family on part of its evidence: its snapshot would say otherwise
                def key = fid.toString()
                if (kinds.size() != expected[key]) {
                    log.warn "Skipping sample QC of ${key}: ${kinds.size()} of its ${expected[key]} inputs were made"
                    return false
                }
                return true
            }
            .map { fid, kinds, payloads ->
                def key = fid.toString()
                def of = { kind -> [kinds, payloads].transpose().findAll { it[0] == kind }.collectMany { it[1] instanceof List ? it[1] : [it[1]] } }
                tuple(key, of('extract'), of('coverage'), of('norm'), qc.snapshots[key], sites, label, classify)
            }

        SAMPLEQC_FAMILY(family_inputs)
        new_family_tables = SAMPLEQC_FAMILY.out.tables
    } else {
        new_family_tables = channel.empty()
    }

    // =====================================================================
    // Step 3: The cohort
    // =====================================================================
    if (qc.need_cohort) {
        def all_samples = qc.extract_members.values().flatten() as Set
        def all_families = qc.extract_members.keySet()

        // Concat, then a single collect, so the cohort waits for this run's own extracts and
        // families. Never publish cohort tables built from part of the cohort.
        cohort_input = new_extracts
            .map { barcode, extract -> tuple('extract', barcode, [extract]) }
            .concat(channel.fromList((all_samples - qc.needed).collect { barcode ->
                tuple('extract', barcode, [extract_file(barcode)]) }))
            .concat(new_family_tables
                .map { fid, samples, pairs, flags, _inputs -> tuple('family', fid, [samples, pairs, flags]) })
            .concat(channel.fromList((all_families - qc.need_family).collect { fid ->
                tuple('family', fid, published_tables(fid)) }))
            .toList()
            .filter { rows ->
                def have = { kind -> rows.findAll { it[0] == kind }.collect { it[1] } as Set }
                def missing = (all_samples - have('extract')) + (all_families - have('family'))
                if (missing) {
                    log.warn "Skipping the cohort's sample QC: nothing for ${missing.sort().take(20).join(', ')}${missing.size() > 20 ? ' and others' : ''}"
                    return false
                }
                return !rows.isEmpty()
            }
            .map { rows ->
                tuple(params.cohort_name,
                      rows.findAll { it[0] == 'extract' }.collectMany { it[2] },
                      rows.findAll { it[0] == 'family' }.collectMany { it[2] },
                      qc.cohort_snapshot, sites, label, classify)
            }

        SAMPLEQC_COHORT(cohort_input)
        cohort_tables_output = SAMPLEQC_COHORT.out.tables
    } else {
        cohort_tables_output = channel.empty()
    }

    emit:
    extracts = new_extracts
    family_tables = new_family_tables
    cohort_tables = cohort_tables_output
}
