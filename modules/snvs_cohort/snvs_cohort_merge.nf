#!/usr/bin/env nextflow

/*
 * SNVs cohort merge module
 * Merges all family common_gt.bcf files into a cohort-level BCF
 */

nextflow.enable.dsl=2

process SNVS_COHORT_MERGE {
  /*
  Merge all family common_gt.bcf files into a cohort-level BCF file

  Parameters
  ----------
  bcf_files : list
    List of tuples [fid, bcf, csi] for all family BCF files
  cohort_name : val
    Name of the cohort

  Returns
  -------
  Tuple of cohort name, merged BCF file, and its index
  */

  tag "$cohort_name"

  publishDir "${params.data}/cohorts/${cohort_name}/vcfs",
    mode: 'copy',
    pattern: "${cohort_name}.common_gt.bcf*"

  container 'docker://staphb/bcftools:latest'
  
  label 'cohort_merge'

  input:
  val cohort_name
  // Staged into a subdirectory so an input can never collide with the output file name.
  // Without this, a cohort named after its only family (cohort "X" containing family "X") stages
  // X.common_gt.bcf alongside an output also called X.common_gt.bcf: bcftools follows the staged
  // symlink and truncates the published family BCF, emitting a header-only cohort file and exit 0.
  path bcf_files, stageAs: 'input_bcfs/*'
  path csi_files, stageAs: 'input_bcfs/*'

  output:
  tuple val(cohort_name), path("${output_bcf}"), path("${output_bcf}.csi"), emit: cohort_bcf

  script:
  output_bcf = "${cohort_name}.common_gt.bcf"

  """
  set -euo pipefail

  # Create temporary file list for bcftools merge (inputs live in their own directory)
  ls input_bcfs/*.bcf > bcf_file_list.txt

  if [ ! -s bcf_file_list.txt ]; then
    echo "ERROR: no family BCF files were staged for cohort ${cohort_name}" >&2
    exit 1
  fi

  # Record counts come from the index, so this stays cheap on large BCFs
  # stderr is suppressed: "cannot determine contig names given the .csi index alone" is benign here
  n_in=0
  while read -r f; do
    n_in=\$(( n_in + \$(bcftools index -n "\$f" 2>/dev/null) ))
  done < bcf_file_list.txt

  # Merge all family BCF files into cohort BCF
  bcftools merge \\
      --threads ${task.cpus} \\
      --merge none \\
      --force-single \\
      --missing-to-ref \\
      --file-list bcf_file_list.txt \\
      --output-type b \\
      --output ${output_bcf}

  # Index the merged BCF
  bcftools index ${output_bcf}

  # bcftools can exit 0 having written only a header, so verify records survived the merge
  n_out=\$(bcftools index -n ${output_bcf} 2>/dev/null)
  echo "Cohort ${cohort_name}: \${n_in} input records across \$(wc -l < bcf_file_list.txt | tr -d '[:space:]') families -> \${n_out} merged records"
  if [ "\${n_in}" -gt 0 ] && [ "\${n_out}" -eq 0 ]; then
    echo "ERROR: cohort merge produced 0 records from \${n_in} input records" >&2
    exit 1
  fi
  """

  stub:
  output_bcf = "${cohort_name}.common_gt.bcf"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  """
}