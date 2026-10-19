#!/usr/bin/env nextflow

/*
 * CNV coverage: cohort module
 * How well each sample of the cohort matches a control panel
 */

nextflow.enable.dsl=2

process CNV_COVERAGE_FIT {
  /*
  Fit every sample of the cohort against a control panel's autosomal median

  One row per (family_id, sample) of the pedigree that has a coverage file: its pedigree
  and coverage sex, its chrX and chrY ratios, and fit_mad, the MAD of
  log2(sample / panel median) over the autosomal bins where the panel median is at least 1.
  genetics-viz flags the samples that match the panel poorly from it. The file depends on
  every sample, so it is rebuilt whenever one is (re)computed.

  Parameters
  ----------
  cohort_name : val
    Name of the cohort
  files : path
    Every sample's bigWig and JSON
  pedigree : path
    The cohort pedigree
  panel_dir : path
    A panel version directory, as cnv_panel writes it
  tag : val
    cov_b<bin>q<mapq>, which the panel's must equal
  panel_version : val
    The panel's panel_version, part of the output file name
  script : path
    cnv_coverage, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of cohort name and the fit table
  */

  tag "${cohort_name}"

  publishDir "${params.data}/cohorts/${cohort_name}/coverage",
    mode: 'copy',
    pattern: "${cohort_name}.${tag}.panel_v${panel_version}.fit.tsv"

  label 'cnv_coverage_fit'

  input:
  // In their own directory, so that a sample named after the cohort cannot collide with an output
  tuple val(cohort_name), path(files, stageAs: 'inputs/*'), path(pedigree), path(panel_dir), val(tag), val(panel_version), path(script)

  output:
  tuple val(cohort_name), path("${cohort_name}.${tag}.panel_v${panel_version}.fit.tsv"), emit: fit

  script:
  """
  set -euo pipefail

  python3 ${script} fit \\
      --pedigree ${pedigree} \\
      --input-dir inputs \\
      --tag ${tag} \\
      --panel-dir ${panel_dir} \\
      --threads ${task.cpus} \\
      --out ${cohort_name}.${tag}.panel_v${panel_version}.fit.tsv
  """

  stub:
  """
  printf 'family_id\\tsample\\tpedigree_sex\\tcoverage_sex\\tfit_mad\\tchrX_ratio\\tchrY_ratio\\n' > ${cohort_name}.${tag}.panel_v${panel_version}.fit.tsv
  """
}
