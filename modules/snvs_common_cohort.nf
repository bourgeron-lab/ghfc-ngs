#!/usr/bin/env nextflow

/*
 * SNVs common cohort merge module
 * Merges all family common_gt.bcf files into a cohort-level BCF
 */

nextflow.enable.dsl=2

process snvs_common_cohort {
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
  path bcf_files
  path csi_files

  output:
  tuple val(cohort_name), path("${output_bcf}"), path("${output_bcf}.csi"), emit: cohort_bcf

  script:
  output_bcf = "${cohort_name}.common_gt.bcf"

  """
  # Create temporary file list for bcftools merge
  ls *.bcf > bcf_file_list.txt
  
  # Merge all family BCF files into cohort BCF
  bcftools merge \\
      --threads ${task.cpus} \\
      --merge none \\
      --missing-to-ref \\
      --file-list bcf_file_list.txt \\
      --output-type b \\
      --output ${output_bcf}

  # Index the merged BCF
  bcftools index ${output_bcf}
  """

  stub:
  output_bcf = "${cohort_name}.common_gt.bcf"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  """
}