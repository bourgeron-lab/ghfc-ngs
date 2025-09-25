#!/usr/bin/env nextflow

/*
 * SNVs common variants filtering module
 * Filters common variants and keeps only GT format field
 */

nextflow.enable.dsl=2

process snvs_common_filters {
  /*
  Filter common variants by removing variants with non-empty INFO/genomes_filters
  and keeping only GT field in FORMAT

  Parameters
  ----------
  fid : str
    Family ID
  vcf : path
    Common variants VCF file (FID.common.vcf.gz)
  vcf_index : path
    VCF index file (.tbi or .csi)

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
  tuple val(fid), path(vcf), path(vcf_index)

  output:
  tuple val(fid), path("${output_bcf}"), path("${output_bcf}.csi"), emit: filtered_common_bcfs

  script:
  output_bcf = "${fid}.common_gt.bcf"

  """
  # Filter common variants: exclude variants with non-empty genomes_filters and keep only GT format field
  bcftools view \\
      --threads ${task.cpus} \\
      -i 'INFO/genomes_filters="." | INFO/genomes_filters=""' \\
      -x '^FORMAT/GT' \\
      -O b \\
      -o ${output_bcf} \\
      ${vcf}

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