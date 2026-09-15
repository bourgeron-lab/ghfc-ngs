#!/usr/bin/env nextflow

/*
 * Ancestry: principal component projection module
 * Projects a family onto the reference panel and assigns an ancestry label
 */

nextflow.enable.dsl=2

process ANCESTRY_PCS {
  /*
  Project a family's panel genotypes onto the reference panel's components

  The projection is per-sample and uses the bundle's own allele frequencies, so a
  family scores identically alone or inside a cohort. The PCs emitted here are the
  ONLY ones the PGS adjustment may consume: projected components shrink toward the
  origin relative to in-sample ones, and the bundle's reference PCs are the
  projected kind, so mixing sources biases the correction.

  Parameters
  ----------
  fid : val
    Family ID
  bcf : path
    Family panel genotype BCF file
  bcf_index : path
    BCF index file (.csi)
  reference : val
    Path to the mounted ancestry-pgs reference bundle
  panel_name : val
    Label identifying this panel/threshold combination in output file names

  Returns
  -------
  Tuple of family ID, PCs table, ancestry label table, and the QC report
  */

  tag "$fid"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/ancestry"
  },
    mode: 'copy',
    pattern: "${fid}.${panel_name}.{pcs.tsv,ancestry.tsv,pcs.qc.json}"

  label 'ancestry_pcs'

  input:
  tuple val(fid), path(bcf), path(bcf_index), val(reference), val(panel_name)

  output:
  tuple val(fid), path("${out_prefix}.pcs.tsv"), emit: pcs
  tuple val(fid), path("${out_prefix}.ancestry.tsv"), emit: ancestry
  tuple val(fid), path("${out_prefix}.pcs.qc.json"), emit: qc

  script:
  out_prefix = "${fid}.${panel_name}"

  """
  set -euo pipefail

  # --npc is deliberately not set: the default takes every component the bundle
  # carries, and the fitted PGS model needs the first ten of them.
  ancestry-pgs pcs \\
      --bcf ${bcf} \\
      --reference ${reference} \\
      --out-prefix ${out_prefix} \\
      --threads ${task.cpus}
  """

  stub:
  out_prefix = "${fid}.${panel_name}"
  """
  printf 'IID\\tPC1\\tPC2\\n${fid}_s1\\t0.001\\t0.002\\n' > ${out_prefix}.pcs.tsv
  printf 'IID\\tregion\\tpopulation\\n${fid}_s1\\tEUR\\tpop_EUR\\n' > ${out_prefix}.ancestry.tsv
  echo '{}' > ${out_prefix}.pcs.qc.json
  """
}
