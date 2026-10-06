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
  set -euo pipefail

  # List all unique DNM TSV files (sorted for consistency)
  # find, not ls: a large cohort expands the glob past the kernel argument limit
  find . -maxdepth 1 -name '*.dnm.tsv' -printf '%f\\n' | sort -u > file_list.txt

  # Refuse to publish an empty cohort file over an existing one
  if [ ! -s file_list.txt ]; then
    echo "ERROR: no family DNM TSV files were staged for cohort ${cohort_name}" >&2
    exit 1
  fi

  # Get the first file to extract the header
  first_file=\$(head -n 1 file_list.txt)

  # Write header from first file
  head -n 1 "\${first_file}" > ${output_tsv}

  # Concatenate all files, skipping their headers (xargs batches the list under the argument limit)
  xargs -d '\\n' tail -q -n +2 < file_list.txt >> ${output_tsv}
  """

  stub:
  output_tsv = "${cohort_name}.${vep_config_name}.dnm.tsv"
  """
  touch ${output_tsv}
  """
}
