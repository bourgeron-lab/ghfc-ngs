#!/usr/bin/env nextflow

/*
 * Wombat cohort merge module
 * Merges all family Wombat TSV files into cohort-level TSV files
 */

nextflow.enable.dsl=2

process MERGE_WOMBAT {
  /*
  Merge all family Wombat TSV files into cohort-level TSV files

  Parameters
  ----------
  cohort_name : val
    Name of the cohort
  wombat_tsv_files : list
    List of Wombat TSV files from all families
  vep_config_name : val
    VEP configuration name for output file naming
  wombat_config_name : val
    Wombat configuration name for output file naming
  output_name : val
    Output file name (from impact section in wombat config)

  Returns
  -------
  Tuple of cohort name, wombat config name, output name, and merged TSV file
  */

  tag "${cohort_name}-${wombat_config_name}-${output_name}"

  publishDir "${params.data}/cohorts/${cohort_name}/wombat",
    mode: 'copy',
    pattern: "${cohort_name}.*.tsv"

  container 'docker://ubuntu:22.04'
  
  label 'merge_wombat'

  input:
  val cohort_name
  path wombat_tsv_files
  val vep_config_name
  val wombat_config_name
  val output_name

  output:
  tuple val(cohort_name), val(wombat_config_name), val(output_name), path("${output_tsv}"), emit: cohort_wombat_tsv

  script:
  output_tsv = "${cohort_name}.rare.${vep_config_name}.${wombat_config_name}.${output_name}.tsv"

  """
  # List all matching TSV files (sorted for consistency)
  # New pattern: {FID}.rare.{vep_config_name}.{wombat_config_name}.tsv
  ls -1 *.${wombat_config_name}.tsv | sort -u > file_list.txt
  
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
  output_tsv = "${cohort_name}.rare.${vep_config_name}.${wombat_config_name}.${output_name}.tsv"
  """
  echo "header" > ${output_tsv}
  """
}
