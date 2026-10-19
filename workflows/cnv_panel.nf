/*
 * CNV panel Workflow
 * Builds a versioned panel of controls - per-bin median and MAD of the controls' normalised
 * coverage, for the autosomes, XX chrX and XY chrX/chrY - from every sample of a controls cohort
 *
 * Built once per cnv_panel_version: an existing version is never rebuilt nor overwritten
 */

// Include modules
include { CNV_PANEL_BUILD } from '../modules/cnv/cnv_panel_build'

workflow CNV_PANEL {

    take:
    new_coverage      // channel: [barcode, bw, json] - coverage computed in this run
    pedigree_file     // string: path to the cohort pedigree
    need_sample       // collection: samples whose coverage this run computes
    existing_samples  // collection: samples whose coverage is already on disk

    main:

    def tag = "cov_b${params.cnv_cov_bin}q${params.cnv_cov_mapq}"
    def script = file("${projectDir}/modules/cnv/scripts/cnv_coverage", checkIfExists: true)

    def expected_samples = [] as Set
    expected_samples.addAll(need_sample)
    expected_samples.addAll(existing_samples)

    // Read from their sharded paths, never globbed
    existing_coverage = channel.fromList(existing_samples.collect { barcode ->
        def prefix = "${Sharding.getSampleDir(params.data, barcode)}/coverage/${barcode}.${tag}"
        tuple(barcode, file("${prefix}.bw"), file("${prefix}.json"))
    })

    // A panel is built once and never rebuilt, so never from part of the cohort: a sample whose
    // coverage failed blocks it, to be built by a later run
    panel_input = new_coverage
        .concat(existing_coverage)
        .unique { barcode, _bw, _json -> barcode }
        .toList()
        .filter { rows ->
            def missing = expected_samples - (rows.collect { barcode, _bw, _json -> barcode } as Set)
            if (missing) {
                log.warn "Skipping the CNV panel build: no coverage for ${missing.sort().take(20).join(', ')}${missing.size() > 20 ? ' and others' : ''}"
                return false
            }
            return !rows.isEmpty()
        }
        .map { rows ->
            tuple(params.cohort_name,
                  rows.collectMany { _barcode, bw, json -> [bw, json] },
                  file(pedigree_file),
                  tag, params.cnv_panel_version.toString(), params.cnv_panel_max_fit_mad, script)
        }

    CNV_PANEL_BUILD(panel_input)

    emit:
    panel = CNV_PANEL_BUILD.out.panel
}
