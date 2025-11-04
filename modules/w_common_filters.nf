#!/usr/bin/env nextflow

/*
 * Wombat: SNVs common variants filtering module
 * Filters common variants and keeps only GT format field
 */

nextflow.enable.dsl=2

process W_COMMON_FILTERS {
  /*
  Filter common variants by removing variants with non-empty INFO/genomes_filters
  and keeping only GT field in FORMAT

  Parameters
  ----------
  fid : str
    Family ID
  bcf : path
    Common variants BCF file (FID.common.bcf)
  bcf_index : path
    BCF index file (.csi)

  Returns
  -------
  Tuple of family ID, filtered BCF file, and its index
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "${fid}.common_gt.bcf*"

  container 'docker://staphb/bcftools:latest'
  
  label 'common_filtering'

  input:
  tuple val(fid), path(bcf), path(bcf_index)

  output:
  tuple val(fid), path("${output_bcf}"), path("${output_bcf}.csi"), emit: filtered_common_bcfs

  script:
  output_bcf = "${fid}.common_gt.bcf"

  """
  # Filter common variants: exclude variants with non-empty genomes_filters and keep only GT format field
  bcftools annotate \\
      --threads ${task.cpus} \\
      -i 'INFO/genomes_filters="." | INFO/genomes_filters=""' \\
      -x '^FORMAT/GT' \\
      -O b \\
      -o ${output_bcf} \\
      ${bcf}

  # Index the filtered BCF
  bcftools index ${output_bcf}
  """

  stub:
  output_bcf = "${fid}.common_gt.bcf"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  """
}
