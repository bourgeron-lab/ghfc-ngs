#!/usr/bin/env nextflow

/*
 * CNV panel: cohort module
 * A versioned panel of controls for the normalised coverage of any cohort
 */

nextflow.enable.dsl=2

process CNV_PANEL_BUILD {
  /*
  Build a control panel from every sample of a controls cohort

  Per bin, the median and the unscaled MAD of the controls' normalised values, for the
  autosomes (every retained control), XX (chrX) and XY (chrX and chrY), after excluding
  the controls that fit the others poorly (fit_mad above max_fit_mad), and keeping out of
  the sex groups those whose pedigree sex and coverage disagree; see the script's
  docstring. The version directory is self-contained, for other cohorts and genetics-viz
  to reference as it is, and is never overwritten: the script refuses an existing one,
  and so does the publication.

  Parameters
  ----------
  cohort_name : val
    Name of the controls cohort
  files : path
    Every sample's bigWig and JSON
  pedigree : path
    The cohort pedigree
  tag : val
    cov_b<bin>q<mapq>
  version : val
    cnv_panel_version; the directory is v<version>
  max_fit_mad : val
    cnv_panel_max_fit_mad
  script : path
    cnv_coverage, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of cohort name and the version directory
  */

  tag "${cohort_name}"

  publishDir "${params.data}/cohorts/${cohort_name}/cnv_panel",
    mode: 'copy',
    overwrite: false

  label 'cnv_panel_build'

  input:
  // In their own directory, so that a sample named after the cohort cannot collide with an output
  tuple val(cohort_name), path(files, stageAs: 'inputs/*'), path(pedigree), val(tag), val(version), val(max_fit_mad), path(script)

  output:
  tuple val(cohort_name), path("v${version}"), emit: panel

  script:
  """
  set -euo pipefail

  python3 ${script} panel \\
      --pedigree ${pedigree} \\
      --input-dir inputs \\
      --tag ${tag} \\
      --cohort ${cohort_name} \\
      --version ${version} \\
      --max-fit-mad ${max_fit_mad} \\
      --threads ${task.cpus} \\
      --out-dir v${version}
  """

  stub:
  """
  mkdir v${version}
  for group in autosomes XX XY; do
    echo stub > v${version}/panel.${tag}.\${group}.median.bw
    echo stub > v${version}/panel.${tag}.\${group}.mad.bw
  done
  printf '{"panel_version": "${version}", "tag": "${tag}"}\\n' > v${version}/panel.json
  printf 'family_id\\tsample\\tpedigree_sex\\tcoverage_sex\\tin_autosomes\\tin_sex_panel\\tautosomal_median_depth\\tchrX_ratio\\tchrY_ratio\\tfit_mad\\texclusion_reason\\n' > v${version}/controls.tsv
  for json in inputs/*.${tag}.json; do
    printf 'FAM\\t%s\\tunknown\\tunknown\\ttrue\\tfalse\\tNA\\tNA\\tNA\\tNA\\t\\n' "\$(basename "\${json}" .${tag}.json)" >> v${version}/controls.tsv
  done
  """
}
