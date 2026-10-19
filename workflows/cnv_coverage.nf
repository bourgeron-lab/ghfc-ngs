/*
 * CNV coverage Workflow
 * Writes each sample's normalised coverage on a fixed bin grid, as a bigWig with a summary JSON,
 * for genetics-viz to show CNVs and SVs against a control panel
 * With cnv_panel_dir, also fits every sample of the cohort against that panel, in one table
 *
 * Skips processing if output files already exist - only runs necessary modules
 */

// Include modules
include { CNV_COVERAGE_SAMPLE } from '../modules/cnv/cnv_coverage_sample'
include { CNV_COVERAGE_FIT } from '../modules/cnv/cnv_coverage_fit'

workflow CNV_COVERAGE {

    take:
    cram_files        // channel: [barcode, cram, crai] - existing and newly aligned CRAMs
    pedigree_file     // string: path to the cohort pedigree
    need_sample       // collection: samples whose coverage must be computed
    existing_samples  // collection: samples whose coverage is already on disk
    need_fit          // boolean: whether the cohort fit table must be (re)built
    panel             // map: [dir, version] of the panel to fit against, or null

    main:

    def tag = "cov_b${params.cnv_cov_bin}q${params.cnv_cov_mapq}"
    def ref = file(params.ref)
    def ref_fai = file("${params.ref}.fai")

    // Staged as a process input rather than referenced through projectDir: the pipeline
    // directory is not mounted inside the task container.
    def script = file("${projectDir}/modules/cnv/scripts/cnv_coverage", checkIfExists: true)

    // =====================================================================
    // Step 1: Each sample's bigWig and JSON
    // =====================================================================
    CNV_COVERAGE_SAMPLE(
        cram_files
            .filter { barcode, _cram, _crai -> barcode in need_sample }
            .unique { barcode, _cram, _crai -> barcode }
            .map { barcode, cram, crai ->
                tuple(barcode, cram, crai, ref, ref_fai, params.cnv_cov_bin, params.cnv_cov_mapq, params.cnv_cov_round, script)
            }
    )
    new_coverage = CNV_COVERAGE_SAMPLE.out.coverage

    // =====================================================================
    // Step 2: The cohort's fit against a panel
    // =====================================================================
    if (need_fit) {
        def expected_samples = [] as Set
        expected_samples.addAll(need_sample)
        expected_samples.addAll(existing_samples)

        // Read from their sharded paths, never globbed
        existing_coverage = channel.fromList(existing_samples.collect { barcode ->
            def prefix = "${Sharding.getSampleDir(params.data, barcode)}/coverage/${barcode}.${tag}"
            tuple(barcode, file("${prefix}.bw"), file("${prefix}.json"))
        })

        // Never publish a table built from part of the cohort: a sample whose coverage failed
        // blocks the fit rather than dropping out of it unnoticed
        fit_input = new_coverage
            .concat(existing_coverage)
            .unique { barcode, _bw, _json -> barcode }
            .toList()
            .filter { rows ->
                def missing = expected_samples - (rows.collect { barcode, _bw, _json -> barcode } as Set)
                if (missing) {
                    log.warn "Skipping the CNV coverage fit: no coverage for ${missing.sort().take(20).join(', ')}${missing.size() > 20 ? ' and others' : ''}"
                    return false
                }
                return !rows.isEmpty()
            }
            .map { rows ->
                tuple(params.cohort_name,
                      rows.collectMany { _barcode, bw, json -> [bw, json] },
                      file(pedigree_file),
                      file(panel.dir),
                      tag, panel.version, script)
            }

        CNV_COVERAGE_FIT(fit_input)
        fit_output = CNV_COVERAGE_FIT.out.fit
    } else {
        fit_output = channel.empty()
    }

    emit:
    coverage = new_coverage   // channel: [barcode, bw, json]
    fit = fit_output
}
