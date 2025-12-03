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
  # List all unique DNM TSV files (sorted for consistency)
  ls -1 *.dnm.tsv | sort -u > file_list.txt
  
  # Get the first file to extract the header
  first_file=\$(head -n 1 file_list.txt)
  
  # Write header from first file
  head -n 1 "\${first_file}" > ${output_tsv}
  
  # Concatenate all files, skipping their headers
  while IFS= read -r file; do
    tail -n +2 "\${file}" >> ${output_tsv}
  done < file_list.txt
  """

  stub:
  output_tsv = "${cohort_name}.${vep_config_name}.dnm.tsv"
  """
  touch ${output_tsv}
  """
}
