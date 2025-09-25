#!/usr/bin/env nextflow

/*
 * SNVs frequency filtering module - separates rare and common variants using gnomAD frequency thresholds
 */

nextflow.enable.dsl=2

process snvs_gnomad_freq_filter {
  /*
  Separate gnomAD-annotated VCF files into rare and common variants based on frequency thresholds

  Parameters
  ----------
  fid : str
    Family ID
  vcf : path
    gnomAD-annotated VCF file
  vcf_index : path
    VCF index file (.tbi or .csi)
  gnomad_filter_field : val
    gnomAD field to filter on (e.g., "AF_genomes")
  gnomad_filter_threshold : val
    Frequency threshold for separation (variants < threshold = rare, >= threshold = common)

  Returns
  -------
  Two tuples: (family ID, rare VCF file, index) and (family ID, common VCF file, index)
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "${fid}.{rare,common}.vcf.gz*"

  container 'docker://staphb/bcftools:latest'
  
  label 'gnomad_filtering'

  input:
  tuple val(fid), path(vcf), path(vcf_index), val(gnomad_filter_field), val(gnomad_filter_threshold)

  output:
  tuple val(fid), path("${rare_vcf}.gz"), path("${rare_vcf}.gz.tbi"), emit: rare_vcfs
  tuple val(fid), path("${common_vcf}.gz"), path("${common_vcf}.gz.tbi"), emit: common_vcfs

  script:
  rare_vcf = "${fid}.rare.vcf"
  common_vcf = "${fid}.common.vcf"

  """
  echo "Filtering VCF for rare variants with ${gnomad_filter_field} < ${gnomad_filter_threshold} or missing"

  # Filter VCF for rare variants (frequency < threshold or missing)
  bcftools view \\
      --threads ${task.cpus} \\
      -i 'INFO/${gnomad_filter_field}<${gnomad_filter_threshold} | INFO/${gnomad_filter_field}="."' \\
      -O z \\
      -o ${rare_vcf}.gz \\
      ${vcf}

  echo "Filtering VCF for common variants with ${gnomad_filter_field} >= ${gnomad_filter_threshold}"

  # Filter VCF for common variants (frequency >= threshold)
  bcftools view \\
      --threads ${task.cpus} \\
      -i 'INFO/${gnomad_filter_field}>=${gnomad_filter_threshold}' \\
      -O z \\
      -o ${common_vcf}.gz \\
      ${vcf}

  echo "Indexing rare and common VCFs"

  # Index both VCFs
  bcftools index --tbi ${rare_vcf}.gz
  bcftools index --tbi ${common_vcf}.gz
  """

  stub:
  rare_vcf = "${fid}.rare.vcf"
  common_vcf = "${fid}.common.vcf"
  """
  touch ${rare_vcf}.gz
  touch ${rare_vcf}.gz.tbi
  touch ${common_vcf}.gz
  touch ${common_vcf}.gz.tbi
  """
}