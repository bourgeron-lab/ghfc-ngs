#!/usr/bin/env nextflow

/*
 * De novo: cohort tables
 * Every family's candidates and QC rows, with recurrence and the cohort's QC range
 */

nextflow.enable.dsl=2

process DENOVO_COHORT {
  /*
  Gather the families' de novo tables for the cohort

  denovo_tables cohort concatenates them and adds, on the candidates, how many other families
  carry the same one (RECURRENT, which makes it not strict: a recurrent site is usually an
  artefact, unless it is a known recurrent pathogenic variant, hence a flag), and on the QC rows
  the strict count recomputed after it, the range (median +- 4 scaled MADs of the rate per
  callable unmasked Gb with 20 children or more, a fixed 15-40 below), and the sampleqc flags of
  the child's trio from the cohort's sampleqc flags table when there is one.

  Parameters
  ----------
  cohort : val
    Cohort name
  tables : path
    Every family's dnm and qc tables
  sampleqc_flags : path
    The cohort's sampleqc flags table, or nothing
  genomeqc_samples : path
    The cohort's genome QC samples table, or nothing: each child's seq_class and CELL_LINE flag
  snapshot : val
    The cohort's inputs snapshot, as text
  label : val
    denovo_label
  tables_script : path
    denovo_tables

  Returns
  -------
  Tuple of cohort name, its dnm and qc tables, and the snapshot
  */

  tag "${cohort}"

  publishDir "${params.data}/cohorts/${cohort}/denovo",
    mode: 'copy',
    pattern: "${cohort}.${label}.{dnm,qc,inputs}.tsv"

  label 'denovo_cohort'

  input:
  tuple val(cohort), path(tables, stageAs: 'families/*'), path(sampleqc_flags, stageAs: 'sampleqc/*'),
        path(genomeqc_samples, stageAs: 'genomeqc/*'), val(snapshot), val(label), path(tables_script)

  output:
  tuple val(cohort), path("${prefix}.dnm.tsv"), path("${prefix}.qc.tsv"), path("${prefix}.inputs.tsv"), emit: tables

  script:
  prefix = "${cohort}.${label}"
  def flags = sampleqc_flags instanceof List ? sampleqc_flags : (sampleqc_flags ? [sampleqc_flags] : [])
  def gq = genomeqc_samples instanceof List ? genomeqc_samples : (genomeqc_samples ? [genomeqc_samples] : [])
  """
  set -euo pipefail

  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}

  python3 ${tables_script} cohort \\
    --dnm families/*.dnm.tsv \\
    --qc families/*.qc.tsv \\
    ${flags ? "--sampleqc-flags ${flags[0]}" : ''} \\
    ${gq ? "--genomeqc-samples ${gq[0]}" : ''} \\
    --out-prefix ${prefix}
  """

  stub:
  prefix = "${cohort}.${label}"
  """
  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}
  printf 'family_id\\tchild\\tchrom\\tpos\\tref\\talt\\tflags\\tstrict\\trecurrence_families\\n' > ${prefix}.dnm.tsv
  printf 'family_id\\tchild\\tn_candidates\\tn_strict\\tqc\\n' > ${prefix}.qc.tsv
  """
}
