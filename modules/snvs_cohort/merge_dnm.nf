#!/usr/bin/env nextflow

/*
 * SNVs cohort DNM merge module
 * Concatenates all family DNM TSV files into a cohort-level TSV
 */

nextflow.enable.dsl=2

process MERGE_DNM {
  /*
  Merge all family DNM TSV files into a cohort-level TSV file

  Parameters
  ----------
  cohort_name : val
    Name of the cohort
  dnm_tsv_files : list
    List of DNM TSV files from all families
  vep_config_name : val
    VEP configuration name for output file naming

  Returns
  -------
  Tuple of cohort name and merged DNM TSV file
  */

  tag "$cohort_name"

  publishDir "${params.data}/cohorts/${cohort_name}/vcfs",
    mode: 'copy',
    pattern: "${cohort_name}.*.dnm.tsv"

  container 'docker://ubuntu:22.04'
  
  label 'merge_dnm'

  input:
  val cohort_name
  path dnm_tsv_files
  val vep_config_name

  output:
  tuple val(cohort_name), path("${output_tsv}"), emit: cohort_dnm_tsv

  script:
  output_tsv = "${cohort_name}.${vep_config_name}.dnm.tsv"

  """
  # Get the first file to extract the header
  first_file=\$(ls *.dnm.tsv | head -n 1)
  
  # Write header from first file
  head -n 1 "\${first_file}" > ${output_tsv}
  
  # Concatenate all files, skipping their headers
  for file in *.dnm.tsv; do
    tail -n +2 "\${file}"
  done >> ${output_tsv}
  """

  stub:
  output_tsv = "${cohort_name}.${vep_config_name}.dnm.tsv"
  """
  touch ${output_tsv}
  """
}
