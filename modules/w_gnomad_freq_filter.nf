#!/usr/bin/env nextflow

/*
 * Wombat: SNVs frequency filtering module - separates rare and common variants using gnomAD frequency thresholds
 */

nextflow.enable.dsl=2

process W_GNOMAD_FREQ_FILTER {
  /*
  Separate gnomAD-annotated BCF files into rare and common variants based on frequency thresholds

  Parameters
  ----------
  fid : str
    Family ID
  bcf : path
    gnomAD-annotated BCF file
  bcf_index : path
    BCF index file (.csi)
  gnomad_filter_field : val
    gnomAD field to filter on (e.g., "AF_genomes")
  gnomad_filter_threshold : val
    Frequency threshold for separation (variants < threshold = rare, >= threshold = common)

  Returns
  -------
  Two tuples: (family ID, rare BCF file, index) and (family ID, common BCF file, index)
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "${fid}.{rare,common}.bcf*"

  container 'docker://staphb/bcftools:latest'
  
  label 'gnomad_filtering'

  input:
  tuple val(fid), path(bcf), path(bcf_index), val(gnomad_filter_field), val(gnomad_filter_threshold)

  output:
  tuple val(fid), path("${rare_bcf}"), path("${rare_bcf}.csi"), emit: rare_bcfs
  tuple val(fid), path("${common_bcf}"), path("${common_bcf}.csi"), emit: common_bcfs

  script:
  rare_bcf = "${fid}.rare.bcf"
  common_bcf = "${fid}.common.bcf"

  """
  echo "Filtering BCF for rare variants with ${gnomad_filter_field} < ${gnomad_filter_threshold} or missing"

  # Filter BCF for rare variants (frequency < threshold or missing)
  bcftools view \\
      --threads ${task.cpus} \\
      -i 'INFO/${gnomad_filter_field}<${gnomad_filter_threshold} | INFO/${gnomad_filter_field}="."' \\
      -O b \\
      -o ${rare_bcf} \\
      ${bcf}

  echo "Filtering BCF for common variants with ${gnomad_filter_field} >= ${gnomad_filter_threshold}"

  # Filter BCF for common variants (frequency >= threshold)
  bcftools view \\
      --threads ${task.cpus} \\
      -i 'INFO/${gnomad_filter_field}>=${gnomad_filter_threshold}' \\
      -O b \\
      -o ${common_bcf} \\
      ${bcf}

  echo "Indexing rare and common BCFs"

  # Index both BCFs
  bcftools index ${rare_bcf}
  bcftools index ${common_bcf}
  """

  stub:
  rare_bcf = "${fid}.rare.bcf"
  common_bcf = "${fid}.common.bcf"
  """
  touch ${rare_bcf}
  touch ${rare_bcf}.csi
  touch ${common_bcf}
  touch ${common_bcf}.csi
  """
}
