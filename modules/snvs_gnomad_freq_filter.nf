#!/usr/bin/env nextflow

/*
 * SNVs frequency filtering module using gnomAD frequency thresholds
 */

nextflow.enable.dsl=2

process snvs_gnomad_freq_filter {
  /*
  Filter gnomAD-annotated VCF files based on frequency thresholds

  Parameters
  ----------
  fid : str
    Family ID
  vcf : path
    gnomAD-annotated VCF file
  vcf_index : path
    VCF index file (.tbi or .csi)
  gnomad_filter_field : val
    gnomAD field to filter on (e.g., "INFO/AF_genomes")
  gnomad_filter_threshold : val
    Frequency threshold for filtering

  Returns
  -------
  Tuple of family ID, filtered VCF file, and its index
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "${fid}.gnomad_filtered.vcf.gz*"

  container 'docker://staphb/bcftools:latest'
  
  label 'gnomad_filtering'

  input:
  tuple val(fid), path(vcf), path(vcf_index), val(gnomad_filter_field), val(gnomad_filter_threshold)

  output:
  tuple val(fid), path("${output_vcf}.gz"), path("${output_vcf}.gz.tbi"), emit: filtered_vcfs

  script:
  output_vcf = "${fid}.gnomad_filtered.vcf"

  """
  # Filter VCF based on gnomAD frequency thresholds
  bcftools view \\
      --threads ${task.cpus} \\
      -i 'INFO/${gnomad_filter_field}<${gnomad_filter_threshold} | INFO/${gnomad_filter_field}="."' \\
      -O z \\
      -o ${output_vcf}.gz \\
      ${vcf}

  # Index the filtered VCF
  bcftools index --tbi ${output_vcf}.gz
  """

  stub:
  output_vcf = "${fid}.gnomad_filtered.vcf"
  """
  touch ${output_vcf}.gz
  touch ${output_vcf}.gz.tbi
  """
}