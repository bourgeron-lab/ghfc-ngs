#!/usr/bin/env nextflow

/*
 * Ancestry: admixture projection module
 * Projects a family onto the reference panel's admixture components
 */

nextflow.enable.dsl=2

process ANCESTRY_ADMIXTURE {
  /*
  Project a family's panel genotypes onto the trained admixture components

  The trained model normalises per sample, so the proportions do not depend on how
  many samples are submitted. It also treats a missing genotype as homozygous
  reference, which biases the proportions instead of erroring, and so the command
  refuses below 90% per-sample coverage of the panel's SNPs. That refusal is the
  reason this step reads gVCF-derived panel genotypes rather than the family's
  common_gt.bcf, which cannot reach that bar.

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
  Tuple of family ID, admixture proportions table, and the QC report
  */

  tag "$fid"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/ancestry"
  },
    mode: 'copy',
    pattern: "${fid}.${panel_name}.{Q.tsv,admixture.qc.json}"

  label 'ancestry_admixture'

  input:
  tuple val(fid), path(bcf), path(bcf_index), val(reference), val(panel_name)

  output:
  tuple val(fid), path("${out_prefix}.Q.tsv"), emit: admixture
  tuple val(fid), path("${out_prefix}.admixture.qc.json"), emit: qc

  script:
  out_prefix = "${fid}.${panel_name}"
  min_coverage_arg = params.ancestry_min_coverage ? "--min-coverage ${params.ancestry_min_coverage}" : ""

  """
  set -euo pipefail

  ancestry-pgs admixture \\
      --bcf ${bcf} \\
      --reference ${reference} \\
      --out-prefix ${out_prefix} \\
      --threads ${task.cpus} \\
      ${min_coverage_arg}
  """

  stub:
  out_prefix = "${fid}.${panel_name}"
  """
  printf 'IID\\tsnp_coverage\\tAFR\\tEUR\\n${fid}_s1\\t0.99\\t0.1\\t0.9\\n' > ${out_prefix}.Q.tsv
  echo '{}' > ${out_prefix}.admixture.qc.json
  """
}
