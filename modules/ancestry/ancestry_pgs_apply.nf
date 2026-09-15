#!/usr/bin/env nextflow

/*
 * Ancestry: polygenic score ancestry adjustment module
 * Applies the fitted z-score model to a family's raw scores
 */

nextflow.enable.dsl=2

process ANCESTRY_PGS_APPLY {
  /*
  Remove each score's ancestry-expected mean, and z-score it against the panel

  pgs-adjusted and pgs-zscore are the two stages of one computation, so they run
  together here rather than as separate processes. Both read the model fitted for
  this catalog and both verify the catalog's checksum against the model's, so a
  catalog that has changed since the fit is refused rather than silently adjusted
  with a stale model.

  Parameters
  ----------
  fid : val
    Family ID
  pgs_raw : path
    Raw scores table from ANCESTRY_PGS_RAW
  pcs : path
    PCs table from ANCESTRY_PCS
  reference : val
    Path to the mounted ancestry-pgs reference bundle
  catalog : val
    Path to the PGS weight catalog the model was fitted from
  panel_name : val
    Label identifying this panel/threshold combination in output file names

  Returns
  -------
  Tuple of family ID, adjusted scores, z-scores, and the two QC reports
  */

  tag "$fid"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/ancestry"
  },
    mode: 'copy',
    pattern: "${fid}.${panel_name}.{pgs_adjusted.tsv,pgs_zscore.tsv,pgs-adjusted.qc.json,pgs-zscore.qc.json}"

  label 'ancestry_pgs_apply'

  input:
  tuple val(fid), path(pgs_raw), path(pcs), val(reference), val(catalog), val(panel_name)

  output:
  tuple val(fid), path("${out_prefix}.pgs_adjusted.tsv"), emit: pgs_adjusted
  tuple val(fid), path("${out_prefix}.pgs_zscore.tsv"), emit: pgs_zscore
  tuple val(fid), path("${out_prefix}.pgs-adjusted.qc.json"), path("${out_prefix}.pgs-zscore.qc.json"), emit: qc

  script:
  out_prefix = "${fid}.${panel_name}"
  model_arg = params.ancestry_model ? "--model ${params.ancestry_model}" : ""

  """
  set -euo pipefail

  ancestry-pgs pgs-adjusted \\
      --reference ${reference} \\
      --pgs ${pgs_raw} \\
      --pcs ${pcs} \\
      --catalog ${catalog} \\
      --out-prefix ${out_prefix} \\
      --threads ${task.cpus} \\
      ${model_arg}

  ancestry-pgs pgs-zscore \\
      --reference ${reference} \\
      --pgs ${pgs_raw} \\
      --pcs ${pcs} \\
      --catalog ${catalog} \\
      --out-prefix ${out_prefix} \\
      --threads ${task.cpus} \\
      ${model_arg}
  """

  stub:
  out_prefix = "${fid}.${panel_name}"
  """
  printf 'IID\\ttrait_a\\n${fid}_s1\\t-0.5\\n' > ${out_prefix}.pgs_adjusted.tsv
  printf 'IID\\ttrait_a\\n${fid}_s1\\t-1.2\\n' > ${out_prefix}.pgs_zscore.tsv
  echo '{}' > ${out_prefix}.pgs-adjusted.qc.json
  echo '{}' > ${out_prefix}.pgs-zscore.qc.json
  """
}
