/*
 * MELT Workflow
 * Calls ALU, LINE1 and SVA mobile element insertions in each sample's CRAM with MELT Single
 * Merges each family's calls into one VCF, and tabulates each child's insertions against its parents' calls
 * Chains the cohort's insertions into loci with their carrier frequency, and ranks every PASS
 * insertion against the cohort, gnomAD-SV, the genes and the NDD gene set
 *
 * Skips processing if output files already exist - only runs necessary modules
 */

// Include modules
include { MELT_CALL } from '../modules/melt/melt_call'
include { MELT_FAMILY } from '../modules/melt/melt_family'
include { MELT_COHORT } from '../modules/melt/melt_cohort'

workflow MELT {

    take:
    cram_files         // channel: [barcode, cram, crai] - existing and newly aligned CRAMs
    family_members     // map: [barcode: fid]
    pedigree_file      // string: path to the cohort pedigree
    need_call          // collection: samples that need a MELT call
    existing_calls     // collection: samples whose call is already on disk
    need_family        // collection: families whose VCF must be (re)built
    need_cohort_merge  // boolean: whether the cohort files must be (re)built
    genes_bed          // string: gencode gene BED (chrom start end symbol)

    main:

    def label = params.melt_label
    def ref = file(params.ref)
    def ref_fai = file("${params.ref}.fai")
    def melt_dir = file(params.melt_dir)

    // Staged as a process input rather than referenced through projectDir: the pipeline
    // directory is not mounted inside the task container.
    def cohort_score = file("${projectDir}/modules/melt/scripts/melt_cohort_score", checkIfExists: true)
    def transmission_script = file("${projectDir}/modules/melt/scripts/melt_family_transmission", checkIfExists: true)

    // A family's VCF holds every member that has, or gets in this run, a call; a member with
    // neither a call nor a CRAM is left out of it, as the plan already warned
    def expected_members = [:]
    family_members.each { barcode, fid ->
        if (fid in need_family && (barcode in need_call || barcode in existing_calls)) {
            expected_members.get(fid, []) << barcode
        }
    }

    // The cohort files hold every sample that has a call
    def expected_samples = [] as Set
    expected_samples.addAll(need_call)
    expected_samples.addAll(existing_calls)

    // =====================================================================
    // Step 1: Call each sample
    // =====================================================================
    MELT_CALL(
        cram_files
            .filter { barcode, _cram, _crai -> barcode in need_call }
            .unique { barcode, _cram, _crai -> barcode }
            .map { barcode, cram, crai -> tuple(barcode, cram, crai, ref, ref_fai, melt_dir, label) }
    )
    new_calls = MELT_CALL.out.vcfs
        .map { barcode, alu, line1, sva -> tuple(barcode, [alu, line1, sva]) }

    // =====================================================================
    // Step 2: Merge each family
    // =====================================================================
    if (need_family) {
        // Members called in an earlier run
        existing_member_vcfs = channel.fromList(
            expected_members.collectMany { _fid, barcodes ->
                barcodes.findAll { barcode -> barcode in existing_calls }.collect { barcode ->
                    tuple(barcode, meltPublishedVcfs(barcode, label))
                }
            })

        // Keyed with the family's expected size, so a family starts as soon as its last member
        // is called rather than when every sample of the cohort is. remainder keeps a family
        // that never fills up, so that the check below can name it instead of it vanishing.
        family_inputs = new_calls
            .mix(existing_member_vcfs)
            .filter { barcode, _vcfs -> family_members[barcode] in need_family }
            .map { barcode, vcfs ->
                def fid = family_members[barcode]
                tuple(groupKey(fid, expected_members[fid].size()), barcode, vcfs)
            }
            .groupTuple(by: 0, remainder: true)
            .filter { fid, barcodes, _vcfs ->
                // Never publish a family VCF built from part of the family
                def missing = expected_members[fid.toString()] - barcodes
                if (missing) {
                    log.warn "Skipping MELT family merge of ${fid}: no call for ${missing.join(', ')}"
                    return false
                }
                return true
            }
            .map { fid, _barcodes, vcfs ->
                tuple(fid.toString(), vcfs.flatten(), file(pedigree_file), label, transmission_script)
            }

        MELT_FAMILY(family_inputs)
        new_family_vcfs = MELT_FAMILY.out.vcf
        new_family_transmission = MELT_FAMILY.out.transmission
    } else {
        new_family_vcfs = channel.empty()
        new_family_transmission = channel.empty()
    }

    // =====================================================================
    // Step 3: Score the cohort
    // =====================================================================
    if (need_cohort_merge) {
        // Samples not called in this run contribute the calls they already have
        existing_sample_vcfs = channel.fromList(
            existing_calls.collect { barcode -> tuple(barcode, meltPublishedVcfs(barcode, label)) })

        // Concat, then a single collect, so the scoring waits for this run's own calls. Never
        // publish cohort files built from part of the cohort: a sample whose call failed
        // blocks the cohort rather than dropping out of its frequencies unnoticed.
        cohort_input = new_calls
            .concat(existing_sample_vcfs)
            .unique { barcode, _vcfs -> barcode }
            .toList()
            .filter { rows ->
                def missing = expected_samples - (rows.collect { barcode, _vcfs -> barcode } as Set)
                if (missing) {
                    log.warn "Skipping MELT cohort scoring: no call for ${missing.sort().take(20).join(', ')}${missing.size() > 20 ? ' and others' : ''}"
                    return false
                }
                return !rows.isEmpty()
            }
            .map { rows ->
                tuple(params.cohort_name,
                      rows.collectMany { _barcode, vcfs -> vcfs },
                      file(params.melt_gnomad_sv),
                      file(params.melt_ndd_genes),
                      params.melt_ndd_columns,
                      file(genes_bed),
                      label, cohort_score)
            }

        MELT_COHORT(cohort_input)
        cohort_tables_output = MELT_COHORT.out.tables
    } else {
        cohort_tables_output = channel.empty()
    }

    emit:
    sample_vcfs = MELT_CALL.out.vcfs
    family_vcfs = new_family_vcfs
    family_transmission = new_family_transmission
    cohort_tables = cohort_tables_output
}

// A sample's three VCFs, read from its sharded path, never globbed. A function rather than a
// closure in the workflow, which the strict syntax parser does not let be called by name.
def meltPublishedVcfs(barcode, label) {
    def prefix = "${Sharding.getSampleDir(params.data, barcode)}/mei/melt/${barcode}.${label}"
    return ['ALU', 'LINE1', 'SVA'].collect { cls -> file("${prefix}.${cls}.vcf") }
}
