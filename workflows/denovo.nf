/*
 * De novo Workflow
 * Gathers the evidence on each family's de novo candidates - the Wombat de novo table of
 * denovo_wombat_config - from its norm.bcf, the members' CRAMs and coverage bedgraphs (each on
 * disk or made in this run) and the mask bundle; phases them with Unfazed; writes the family's
 * tables, then the cohort's.
 *
 * Rebuilds a family when its inputs snapshot - its pedigree rows, which members have a CRAM and
 * a bedgraph, the Wombat table and mask bundle used - is not the one on disk, when its Wombat
 * table or BCF is remade in this run or is newer than its tables, and the cohort when any family
 * is rebuilt or its snapshot changes
 */

// Include modules
include { DENOVO_EVIDENCE } from '../modules/denovo/denovo_evidence'
include { DENOVO_ORIGIN } from '../modules/denovo/denovo_origin'
include { DENOVO_FAMILY } from '../modules/denovo/denovo_family'
include { DENOVO_COHORT } from '../modules/denovo/denovo_cohort'

workflow DENOVO {

    take:
    norm_bcfs         // channel: [fid, bcf, csi] - existing and newly called family BCFs
    new_wombat        // channel: [fid, config_name, tsv] - Wombat tables made in this run
    new_crams         // channel: [barcode, cram, crai] - CRAMs aligned in this run
    new_bedgraphs     // channel: [barcode, bedgraph, tbi] - coverage bedgraphs made in this run
    new_sampleqc      // channel: flags.tsv - the cohort's sampleqc flags table, when this run writes it
    family_members    // map: [barcode: fid]
    denovo            // map: the plan's denovo entry; see createAnalysisPlan

    main:

    def label = params.denovo_label
    def config_name = params.denovo_wombat_config

    // Staged as a process input rather than referenced through projectDir: the pipeline
    // directory is not mounted inside the task container.
    def tables_script = file("${projectDir}/modules/denovo/scripts/denovo_tables", checkIfExists: true)
    def masks = file(params.denovo_masks.toString(), checkIfExists: true)

    // =====================================================================
    // Step 1: Each family
    // =====================================================================
    if (denovo.need_family) {
        // Every input of a family, tagged with its kind: the Wombat table, the BCF, and one CRAM
        // and one bedgraph per member the snapshot says has one
        def expected = [:]
        denovo.need_family.each { fid ->
            expected[fid] = 2 + denovo.crams[fid].size() + denovo.bedgraphs[fid].size()
        }
        def from_disk = []
        denovo.need_family.each { fid ->
            if (!(fid in denovo.wombat_from_run)) {
                from_disk << tuple(fid, 'table', denovoWombatTable(fid))
            }
            denovo.crams[fid].findAll { barcode -> !(barcode in denovo.crams_from_run) }.each { barcode ->
                from_disk << tuple(fid, 'cram', denovoCram(barcode))
            }
            denovo.bedgraphs[fid].findAll { barcode -> !(barcode in denovo.bedgraphs_from_run) }.each { barcode ->
                from_disk << tuple(fid, 'bedgraph', denovoBedgraph(barcode))
            }
        }

        family_items = channel.fromList(from_disk)
            .mix(new_wombat
                .filter { fid, config, _tsv -> config == config_name && fid in denovo.need_family &&
                          fid in denovo.wombat_from_run }
                .map { fid, _config, tsv -> tuple(fid, 'table', tsv) })
            .mix(norm_bcfs
                .filter { fid, _bcf, _csi -> fid in denovo.need_family }
                .unique { fid, _bcf, _csi -> fid }
                .map { fid, bcf, csi -> tuple(fid, 'norm', [bcf, csi]) })
            .mix(new_crams
                .filter { barcode, _cram, _crai -> barcode in denovo.crams_from_run &&
                          family_members[barcode] in denovo.need_family }
                .map { barcode, cram, crai -> tuple(family_members[barcode], 'cram', [cram, crai]) })
            .mix(new_bedgraphs
                .filter { barcode, _bg, _tbi -> barcode in denovo.bedgraphs_from_run &&
                          family_members[barcode] in denovo.need_family }
                .map { barcode, bg, _tbi -> tuple(family_members[barcode], 'bedgraph', bg) })

        // Keyed with the family's expected count, so a family starts as soon as its last input
        // is made. remainder keeps one that never fills up, so that the check below names it.
        family_inputs = family_items
            .map { fid, kind, payload -> tuple(groupKey(fid, expected[fid]), kind, payload) }
            .groupTuple(by: 0, remainder: true)
            .filter { fid, kinds, _payloads ->
                // Never analyse a family on part of its evidence: its snapshot would say otherwise
                def key = fid.toString()
                if (kinds.size() != expected[key]) {
                    log.warn "Skipping the de novo analysis of ${key}: ${kinds.size()} of its ${expected[key]} inputs were made"
                    return false
                }
                return true
            }
            .map { fid, kinds, payloads ->
                def key = fid.toString()
                tuple(key, denovoFilesOf(kinds, payloads, 'table'), denovoFilesOf(kinds, payloads, 'norm'),
                      denovoFilesOf(kinds, payloads, 'cram'), denovoFilesOf(kinds, payloads, 'bedgraph'))
            }

        evidence_input = family_inputs
            .map { fid, table, norm, crams, bedgraphs ->
                tuple(fid, table, norm, crams, bedgraphs, denovo.snapshots[fid], masks, label, tables_script)
            }
        DENOVO_EVIDENCE(evidence_input)

        origin_input = DENOVO_EVIDENCE.out.evidence
            .map { fid, _evidence, _callable, unfazed, _inputs -> tuple(fid, unfazed) }
            .join(family_inputs.map { fid, _table, _norm, crams, _bedgraphs -> tuple(fid, crams) })
            .map { fid, unfazed, crams -> tuple(fid, unfazed, crams, label, tables_script) }
        DENOVO_ORIGIN(origin_input)

        family_input = DENOVO_EVIDENCE.out.evidence
            .join(DENOVO_ORIGIN.out.origin)
            .map { fid, evidence, callable, unfazed, inputs, origin ->
                tuple(fid, evidence, callable, unfazed, inputs, origin, label, tables_script)
            }
        DENOVO_FAMILY(family_input)
        new_family_tables = DENOVO_FAMILY.out.tables
    } else {
        new_family_tables = channel.empty()
    }

    // =====================================================================
    // Step 2: The cohort
    // =====================================================================
    if (denovo.need_cohort) {
        def all_families = denovo.snapshots.keySet()
        // This run's sampleqc flags table, waited for, or the one on disk, or none. An if, not a
        // ternary: Nextflow 24.10's parser reads `: channel.value(...)` as declaring `channel`
        def disk_flags = denovo.sampleqc_flags ? [file(denovo.sampleqc_flags)] : []
        def flags = new_sampleqc.toList().map { files -> files.take(1) }
        if (!denovo.sampleqc_from_run) {
            flags = channel.value(disk_flags)
        }

        // Concat, then a single collect, so the cohort waits for this run's own families.
        // Never publish cohort tables built from part of the cohort.
        cohort_input = new_family_tables
            .map { fid, dnm, qc, _inputs -> tuple(fid, [dnm, qc]) }
            .concat(channel.fromList((all_families - denovo.need_family).collect { fid ->
                tuple(fid, denovoFamilyTables(fid)) }))
            .toList()
            .filter { rows ->
                def missing = all_families - (rows.collect { row -> row[0] } as Set)
                if (missing) {
                    log.warn "Skipping the cohort's de novo tables: nothing for ${missing.sort().take(20).join(', ')}${missing.size() > 20 ? ' and others' : ''}"
                    return false
                }
                return !rows.isEmpty()
            }
            // Both sides wrapped: combine spreads a list's items
            .map { rows -> [rows] }
            .combine(flags.map { files -> [files] })
            .map { rows, qc_flags ->
                tuple(params.cohort_name, rows.collectMany { row -> row[1] }, qc_flags, denovo.cohort_snapshot,
                      label, tables_script)
            }

        DENOVO_COHORT(cohort_input)
        cohort_tables_output = DENOVO_COHORT.out.tables
    } else {
        cohort_tables_output = channel.empty()
    }

    emit:
    family_tables = new_family_tables
    cohort_tables = cohort_tables_output
}

// Read from their sharded paths, never globbed
def denovoWombatTable(fid) {
    return file("${Sharding.getFamilyDir(params.data, fid)}/wombat/${fid}.rare.${params.vep_config_name}.annotated.${params.denovo_wombat_config}.tsv")
}

def denovoCram(barcode) {
    def cram = file("${Sharding.getSampleDir(params.data, barcode)}/sequences/${barcode}.${params.ref_name}.cram")
    return [cram, file("${cram}.crai")]
}

def denovoBedgraph(barcode) {
    return file("${Sharding.getSampleDir(params.data, barcode)}/sequences/${barcode}.by${params.bin}.bedgraph.gz")
}

// A grouped family's files of one kind, from the kinds and payloads groupTuple collected; the
// BCF's and the CRAMs' payloads are pairs of files
def denovoFilesOf(kinds, payloads, kind) {
    return [kinds, payloads].transpose()
        .findAll { pair -> pair[0] == kind }
        .collectMany { pair -> pair[1] instanceof List ? pair[1] : [pair[1]] }
}

def denovoFamilyTables(fid) {
    def prefix = "${Sharding.getFamilyDir(params.data, fid)}/denovo/${fid}.${params.denovo_label}"
    return ['dnm', 'qc'].collect { kind -> file("${prefix}.${kind}.tsv") }
}
