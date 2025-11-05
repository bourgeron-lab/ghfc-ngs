#!/usr/bin/env nextflow

/*
 * Wombat: SNVs frequency filtering module - separates rare and common variants using gnomAD frequency thresholds
 */

nextflow.enable.dsl=2

process GNOMAD_FREQ_FILTER {
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
    pattern: "${fid}.{rare.vcf.gz*,common.bcf*}"

  container 'docker://staphb/bcftools:latest'
  
  label 'gnomad_filtering'

  input:
  tuple val(fid), path(bcf), path(bcf_index), val(gnomad_filter_field), val(gnomad_filter_threshold)

  output:
  tuple val(fid), path("${rare_vcf}"), path("${rare_vcf}.tbi"), emit: rare_vcfs
  tuple val(fid), path("${common_bcf}"), path("${common_bcf}.csi"), emit: common_bcfs

  script:
  rare_vcf = "${fid}.rare.vcf.gz"
  common_bcf = "${fid}.common.bcf"

  """
  echo "Filtering BCF for rare variants with ${gnomad_filter_field} < ${gnomad_filter_threshold} or missing"

  # Filter BCF for rare variants (frequency < threshold or missing) and output as VCF.gz
  bcftools view \\
      --threads ${task.cpus} \\
      -i 'INFO/${gnomad_filter_field}<${gnomad_filter_threshold} | INFO/${gnomad_filter_field}="."' \\
      -O z \\
      -o ${rare_vcf} \\
      ${bcf}

  echo "Filtering BCF for common variants with ${gnomad_filter_field} >= ${gnomad_filter_threshold}"

  # Filter BCF for common variants (frequency >= threshold)
  bcftools view \\
      --threads ${task.cpus} \\
      -i 'INFO/${gnomad_filter_field}>=${gnomad_filter_threshold}' \\
      -O b \\
      -o ${common_bcf} \\
      ${bcf}

  echo "Indexing rare VCF.gz and common BCF"

  # Index rare VCF.gz with tabix
  tabix -p vcf ${rare_vcf}
  
  # Index common BCF
  bcftools index ${common_bcf}
  """

  stub:
  rare_vcf = "${fid}.rare.vcf.gz"
  common_bcf = "${fid}.common.bcf"
  """
  touch ${rare_vcf}
  touch ${rare_vcf}.tbi
  touch ${common_bcf}
  touch ${common_bcf}.csi
  """
}
