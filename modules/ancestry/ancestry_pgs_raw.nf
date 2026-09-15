#!/usr/bin/env nextflow

/*
 * Ancestry: raw polygenic score module
 * Scores every catalog trait for a family in one pass
 */

nextflow.enable.dsl=2

process ANCESTRY_PGS_RAW {
  /*
  Compute raw polygenic scores for all catalog traits

  Scoring reads the reference bundle's allele frequencies rather than deriving them
  from the samples being scored. That is what makes a family scorable at all:
  PLINK 2 refuses to derive frequencies below 50 samples, and above 50 it would
  quietly use the cohort's own, making the scores incomparable with the reference
  the z-scores are calibrated against.

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
  catalog : val
    Path to the PGS weight catalog (SbayesRC layout)
  panel_name : val
    Label identifying this panel/threshold combination in output file names

  Returns
  -------
  Tuple of family ID, raw scores table, and the QC report
  */

  tag "$fid"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/ancestry"
  },
    mode: 'copy',
    pattern: "${fid}.${panel_name}.{pgs_raw.tsv,pgs-raw.qc.json}"

  label 'ancestry_pgs_raw'

  input:
  tuple val(fid), path(bcf), path(bcf_index), val(reference), val(catalog), val(panel_name)

  output:
  tuple val(fid), path("${out_prefix}.pgs_raw.tsv"), emit: pgs_raw
  tuple val(fid), path("${out_prefix}.pgs-raw.qc.json"), emit: qc

  script:
  out_prefix = "${fid}.${panel_name}"

  """
  set -euo pipefail

  ancestry-pgs pgs-raw \\
      --bcf ${bcf} \\
      --reference ${reference} \\
      --catalog ${catalog} \\
      --out-prefix ${out_prefix} \\
      --threads ${task.cpus}
  """

  stub:
  out_prefix = "${fid}.${panel_name}"
  """
  printf 'IID\\tallele_ct\\ttrait_a\\n${fid}_s1\\t1000\\t0.5\\n' > ${out_prefix}.pgs_raw.tsv
  echo '{}' > ${out_prefix}.pgs-raw.qc.json
  """
}
