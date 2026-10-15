#!/usr/bin/env nextflow

/*
 * MELT: cohort module
 * The cohort frequency of every insertion locus, and every PASS insertion ranked
 */

nextflow.enable.dsl=2

process MELT_COHORT {
  /*
  Score every sample's PASS insertions against the cohort and gnomAD-SV

  The gnomAD-SV mobile element sites are read with bcftools, then melt_cohort_score
  chains every sample's PASS sites into cohort loci with their carrier frequency, and
  ranks every PASS insertion by its composite score; see the script's docstring. Both
  depend on the whole cohort, so the table is rebuilt whenever a sample is (re)called.

  Parameters
  ----------
  cohort_name : val
    Name of the cohort
  vcfs : path
    Every sample's ALU, LINE1 and SVA MELT VCFs
  gnomad_sv : path
    gnomAD-SV sites VCF
  ndd_genes : path
    NDD gene table, with a Symbol column
  ndd_columns : val
    '|'-separated columns of that table which mark a gene as NDD
  genes_bed : path
    Gencode gene BED (chrom start end symbol)
  label : val
    melt_label, identifying the calls in input and output file names
  cohort_score : path
    The scoring script, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of cohort name, the loci BED and the candidates table
  */

  tag "${cohort_name}"

  publishDir "${params.data}/cohorts/${cohort_name}/mei/melt",
    mode: 'copy',
    pattern: "${cohort_name}.${label}{.internal_af.bed,.candidates.tsv.gz}"

  label 'melt_cohort'

  input:
  // In their own directory, so that a sample named after the cohort cannot collide with an output
  tuple val(cohort_name), path(vcfs, stageAs: 'input_vcfs/*'), path(gnomad_sv), path(ndd_genes), val(ndd_columns), path(genes_bed), val(label), path(cohort_score)

  output:
  tuple val(cohort_name), path("${prefix}.internal_af.bed"), path("${prefix}.candidates.tsv.gz"), emit: tables

  script:
  prefix = "${cohort_name}.${label}"

  """
  set -euo pipefail

  # Only the mobile element insertions; the scoring script reads them as written here
  bcftools view -i 'ALT ~ "INS:ME"' ${gnomad_sv} \\
    | bcftools query -f '%CHROM\\t%POS\\t%ALT\\t%INFO/AF\\t%INFO/GRPMAX_AF\\n' > gnomad_mei.tsv
  if [ ! -s gnomad_mei.tsv ]; then
    echo "ERROR: ${gnomad_sv} has no <INS:ME:*> site" >&2
    exit 1
  fi

  python3 ${cohort_score} \\
      --input-dir input_vcfs \\
      --cohort ${cohort_name} \\
      --label ${label} \\
      --genes ${genes_bed} \\
      --ndd ${ndd_genes} \\
      --ndd-columns '${ndd_columns}' \\
      --gnomad gnomad_mei.tsv
  """

  stub:
  prefix = "${cohort_name}.${label}"
  """
  printf 'chr1\\t999\\t1000\\tALU\\t1\\t1\\t1\\n' > ${prefix}.internal_af.bed
  printf 'sample\\tchrom\\tpos\\n' | gzip > ${prefix}.candidates.tsv.gz
  """
}
