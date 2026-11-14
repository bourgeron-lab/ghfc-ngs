#!/usr/bin/env nextflow

/*
 * Genome QC: cohort module
 * Every sample's metrics and variant profile, its seq_class and its flags
 */

nextflow.enable.dsl=2

process GENOMEQC_COHORT {
  /*
  Gather the cohort's sample and family tables, set each sample's seq_class and flag it

  seq_class is the read length class (rl100 under 125 bp, rl150 above), with _pcr for a PCR
  library according to genomeqc_sample_info. Outliers are tested within the class only, against
  its median and scaled MAD, so an older 100 bp genome is never flagged for being one. CELL_LINE comes from the EBV
  genomes per cell (or a dna_source of cell_line); see genomeqc_tables for every flag.

  Parameters
  ----------
  cohort : val
    Cohort name
  samples : path
    The sample tables
  families : path
    The family tables
  pedigree : path
    The cohort's pedigree, for the family of a sample without a family table
  sample_info : path
    genomeqc_sample_info, or nothing
  snapshot : val
    The cohort's inputs snapshot, as text
  label : val
    genomeqc_label
  tables_script : path
    genomeqc_tables

  Returns
  -------
  Tuple of cohort name, its samples and flags tables, and the snapshot
  */

  tag "${cohort}"

  publishDir "${params.data}/cohorts/${cohort}/qc",
    mode: 'copy',
    pattern: "${cohort}.${label}.{samples,flags,inputs}.tsv"

  label 'genomeqc_cohort'

  input:
  tuple val(cohort), path(samples, stageAs: 'samples/*'), path(families, stageAs: 'families/*'),
        path(pedigree, stageAs: 'pedigree/*'), path(sample_info, stageAs: 'info/*'), val(snapshot), val(label),
        path(tables_script)

  output:
  tuple val(cohort), path("${prefix}.samples.tsv"), path("${prefix}.flags.tsv"), path("${prefix}.inputs.tsv"),
        emit: tables

  script:
  prefix = "${cohort}.${label}"
  def info = sample_info instanceof List ? sample_info : (sample_info ? [sample_info] : [])
  """
  set -euo pipefail

  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}

  mkdir -p samples families
  # Directories rather than the paths themselves, which a large cohort would take past ARG_MAX
  python3 ${tables_script} cohort \\
      --samples samples \\
      --families families \\
      --pedigree ${pedigree} \\
      ${info ? "--sample-info ${info[0]}" : ''} \\
      --out-prefix ${prefix}
  """

  stub:
  prefix = "${cohort}.${label}"
  """
  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}
  printf 'sample\\tfamily_id\\tseq_class\\tflags\\n' > ${prefix}.samples.tsv
  printf 'sample\\tfamily_id\\tseq_class\\tflag\\tseverity\\tevidence\\n' > ${prefix}.flags.tsv
  """
}
