/*
 * UPD Workflow
 * Finds each family's uniparental disomy segments from its norm.bcf: with both parents, one, or
 * none (from the roh step's ROH segments), with the cnv_coverage JSONs for the karyotype when
 * there are
 * Gathers the family tables for the cohort, marking the events other families share
 *
 * Rebuilds a family when its inputs snapshot - its pedigree rows and the evidence used - is not
 * the one on disk, or when its BCF or ROH tables are remade in this run, and the cohort when
 * any family is rebuilt or its snapshot changes
 */

// Include modules
include { UPD_FAMILY } from '../modules/upd/upd_family'
include { UPD_COHORT } from '../modules/upd/upd_cohort'

workflow UPD {

    take:
    norm_bcfs         // channel: [fid, bcf, csi] - existing and newly called family BCFs
    new_coverage      // channel: [barcode, bw, json] - coverage the cnv_coverage step makes in this run
    new_roh           // channel: [fid, roh, froh] - ROH tables the roh step makes in this run
    family_members    // map: [barcode: fid]
    upd               // map: the plan's upd entry; see createAnalysisPlan

    main:

    def label = params.upd_label

    // Staged as a process input rather than referenced through projectDir: the pipeline
    // directory is not mounted inside the task container.
    def segments = file("${projectDir}/modules/upd/scripts/upd_segments", checkIfExists: true)

    // =====================================================================
    // Step 1: Each family
    // =====================================================================
    if (upd.need_family) {
        // Every input of a family, tagged with its kind: the BCF, one coverage JSON per member
        // the snapshot says has one, and the ROH tables when it says so
        def expected = [:]
        upd.need_family.each { fid ->
            expected[fid] = 1 + upd.coverage_members[fid].size() + (fid in upd.roh_families ? 1 : 0)
        }
        def from_disk = []
        upd.need_family.each { fid ->
            upd.coverage_members[fid].findAll { barcode -> !(barcode in upd.coverage_from_run) }.each { barcode ->
                from_disk << tuple(fid, 'coverage', updCoverageFile(barcode))
            }
            if (fid in upd.roh_families && !(fid in upd.roh_from_run)) {
                from_disk << tuple(fid, 'roh', updRohFiles(fid))
            }
        }

        family_items = channel.fromList(from_disk)
            .mix(norm_bcfs
                .filter { fid, _bcf, _csi -> fid in upd.need_family }
                .unique { fid, _bcf, _csi -> fid }
                .map { fid, bcf, csi -> tuple(fid, 'norm', [bcf, csi]) })
            .mix(new_coverage
                .filter { barcode, _bw, _json -> barcode in upd.coverage_from_run &&
                          family_members[barcode] in upd.need_family }
                .map { barcode, _bw, json -> tuple(family_members[barcode], 'coverage', json) })
            .mix(new_roh
                .filter { fid, _roh, _froh -> fid in upd.need_family && fid in upd.roh_from_run }
                .map { fid, roh, froh -> tuple(fid, 'roh', [roh, froh]) })

        // Keyed with the family's expected count, so a family starts as soon as its last input
        // is made. remainder keeps one that never fills up, so that the check below names it.
        family_inputs = family_items
            .map { fid, kind, payload -> tuple(groupKey(fid, expected[fid]), kind, payload) }
            .groupTuple(by: 0, remainder: true)
            .filter { fid, kinds, _payloads ->
                // Never analyse a family on part of its evidence: its snapshot would say otherwise
                def key = fid.toString()
                if (kinds.size() != expected[key]) {
                    log.warn "Skipping the UPD analysis of ${key}: ${kinds.size()} of its ${expected[key]} inputs were made"
                    return false
                }
                return true
            }
            .map { fid, kinds, payloads ->
                def key = fid.toString()
                tuple(key, updFilesOf(kinds, payloads, 'norm'), updFilesOf(kinds, payloads, 'coverage'),
                      updFilesOf(kinds, payloads, 'roh'), upd.snapshots[key], params.gnomad_file, upd.mask ?: '',
                      label, segments)
            }

        UPD_FAMILY(family_inputs)
        new_family_tables = UPD_FAMILY.out.tables
    } else {
        new_family_tables = channel.empty()
    }

    // =====================================================================
    // Step 2: The cohort
    // =====================================================================
    if (upd.need_cohort) {
        def all_families = upd.snapshots.keySet()

        // Concat, then a single collect, so the cohort waits for this run's own families.
        // Never publish cohort tables built from part of the cohort.
        cohort_input = new_family_tables
            .map { fid, seg, events, members, _windows, _inputs -> tuple(fid, [seg, events, members]) }
            .concat(channel.fromList((all_families - upd.need_family).collect { fid ->
                tuple(fid, updFamilyTables(fid)) }))
            .toList()
            .filter { rows ->
                def missing = all_families - (rows.collect { row -> row[0] } as Set)
                if (missing) {
                    log.warn "Skipping the cohort's UPD tables: nothing for ${missing.sort().take(20).join(', ')}${missing.size() > 20 ? ' and others' : ''}"
                    return false
                }
                return !rows.isEmpty()
            }
            .map { rows ->
                tuple(params.cohort_name, rows.collectMany { row -> row[1] }, upd.cohort_snapshot, label, segments)
            }

        UPD_COHORT(cohort_input)
        cohort_tables_output = UPD_COHORT.out.tables
    } else {
        cohort_tables_output = channel.empty()
    }

    emit:
    family_tables = new_family_tables
    cohort_tables = cohort_tables_output
}

// Read from their sharded paths, never globbed
def updCoverageFile(barcode) {
    return file("${Sharding.getSampleDir(params.data, barcode)}/coverage/${barcode}.cov_b${params.cnv_cov_bin}q${params.cnv_cov_mapq}.json")
}

def updRohFiles(fid) {
    def prefix = "${Sharding.getFamilyDir(params.data, fid)}/roh/${fid}.${params.ancestry_panel_name}"
    return [file("${prefix}.roh.tsv"), file("${prefix}.froh.tsv")]
}

def updFamilyTables(fid) {
    def prefix = "${Sharding.getFamilyDir(params.data, fid)}/upd/${fid}.${params.upd_label}"
    return ['segments', 'events', 'members'].collect { kind -> file("${prefix}.${kind}.tsv") }
}

// A grouped family's files of one kind, from the kinds and payloads groupTuple collected; the
// BCF's and the ROH tables' payloads are pairs of files
def updFilesOf(kinds, payloads, kind) {
    return [kinds, payloads].transpose()
        .findAll { pair -> pair[0] == kind }
        .collectMany { pair -> pair[1] instanceof List ? pair[1] : [pair[1]] }
}
