/*
 * Genome QC Workflow
 * Measures what kind of DNA and library each sample is, from its CRAM: EBV and mtDNA copies per
 * cell, read length, duplicates, mismatches, insert size and indels in the alignments
 * Profiles each family member's SNVs and indels from the family BCF
 * Gathers both for the cohort, gives each sample its seq_class (read length, PCR) and flags the
 * cell lines and the outliers of their class
 *
 * A sample is measured once, when its table is missing. Rebuilds a family when its inputs
 * snapshot - its pedigree rows and the TR BED - is not the one on disk, or when its BCF is
 * called in this run or is newer than its table, and the cohort when any sample or family is
 * redone or its snapshot changes
 */

// Include modules
include { GENOMEQC_SAMPLE } from '../modules/genomeqc/genomeqc_sample'
include { GENOMEQC_FAMILY } from '../modules/genomeqc/genomeqc_family'
include { GENOMEQC_COHORT } from '../modules/genomeqc/genomeqc_cohort'

workflow GENOMEQC {

    take:
    cram_files        // channel: [barcode, cram, crai] - existing and newly aligned CRAMs
    norm_bcfs         // channel: [fid, bcf, csi] - existing and newly called family BCFs
    gq                // map: the plan's genomeqc entry; see createAnalysisPlan

    main:

    def label = params.genomeqc_label
    def ref = file(params.ref)
    def ref_fai = file("${params.ref}.fai")
    // Staged as a process input rather than referenced through projectDir: the pipeline
    // directory is not mounted inside the task container.
    def tables_script = file("${projectDir}/modules/genomeqc/scripts/genomeqc_tables", checkIfExists: true)
    def tr_bed = params.genomeqc_tr_bed ? [file(params.genomeqc_tr_bed.toString(), checkIfExists: true)] : []
    def sample_info = params.genomeqc_sample_info ? [file(params.genomeqc_sample_info.toString(), checkIfExists: true)] : []

    // =====================================================================
    // Step 1: Each sample
    // =====================================================================
    GENOMEQC_SAMPLE(
        cram_files
            .filter { barcode, _cram, _crai -> barcode in gq.needed }
            .unique { barcode, _cram, _crai -> barcode }
            .map { barcode, cram, crai -> tuple(barcode, cram, crai, ref, ref_fai, label, tables_script) }
    )
    new_samples = GENOMEQC_SAMPLE.out.metrics

    // =====================================================================
    // Step 2: Each family with a BCF
    // =====================================================================
    if (gq.need_family) {
        GENOMEQC_FAMILY(
            norm_bcfs
                .filter { fid, _bcf, _csi -> fid in gq.need_family }
                .unique { fid, _bcf, _csi -> fid }
                .map { fid, bcf, csi -> tuple(fid, [bcf, csi], tr_bed, gq.snapshots[fid], label, tables_script) }
        )
        new_family_tables = GENOMEQC_FAMILY.out.tables
    } else {
        new_family_tables = channel.empty()
    }

    // =====================================================================
    // Step 3: The cohort
    // =====================================================================
    if (gq.need_cohort) {
        def all_samples = (gq.existing + gq.needed) as Set
        def all_families = gq.snapshots.keySet()
        def pedigree = file(RunContext.state.pedigree_file.toString())

        // Concat, then a single collect, so the cohort waits for this run's own samples and
        // families. Never publish cohort tables built from part of the cohort.
        cohort_input = new_samples
            .map { barcode, metrics, _stats -> tuple('sample', barcode, [metrics]) }
            .concat(channel.fromList((all_samples - gq.needed).collect { barcode ->
                tuple('sample', barcode, [genomeqcSampleTable(barcode)]) }))
            .concat(new_family_tables
                .map { fid, variants, _inputs -> tuple('family', fid, [variants]) })
            .concat(channel.fromList((all_families - gq.need_family).collect { fid ->
                tuple('family', fid, [genomeqcFamilyTable(fid)]) }))
            .toList()
            .filter { rows ->
                def missing = (all_samples - genomeqcIdsOf(rows, 'sample')) + (all_families - genomeqcIdsOf(rows, 'family'))
                if (missing) {
                    log.warn "Skipping the cohort's genome QC: nothing for ${missing.sort().take(20).join(', ')}${missing.size() > 20 ? ' and others' : ''}"
                    return false
                }
                return !rows.isEmpty()
            }
            .map { rows ->
                tuple(params.cohort_name,
                      rows.findAll { row -> row[0] == 'sample' }.collectMany { row -> row[2] },
                      rows.findAll { row -> row[0] == 'family' }.collectMany { row -> row[2] },
                      pedigree, sample_info, gq.cohort_snapshot, label, tables_script)
            }

        GENOMEQC_COHORT(cohort_input)
        cohort_tables_output = GENOMEQC_COHORT.out.tables
    } else {
        cohort_tables_output = channel.empty()
    }

    emit:
    sample_tables = new_samples
    family_tables = new_family_tables
    cohort_tables = cohort_tables_output
}

// Read from their sharded paths, never globbed
def genomeqcSampleTable(barcode) {
    return file("${Sharding.getSampleDir(params.data, barcode)}/qc/genomeqc/${barcode}.${params.genomeqc_label}.tsv")
}

def genomeqcFamilyTable(fid) {
    return file("${Sharding.getFamilyDir(params.data, fid)}/qc/${fid}.${params.genomeqc_label}.variants.tsv")
}

// The IDs of one kind among the cohort's [kind, id, files] rows
def genomeqcIdsOf(rows, kind) {
    return rows.findAll { row -> row[0] == kind }.collect { row -> row[1] } as Set
}
