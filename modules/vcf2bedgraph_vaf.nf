#!/usr/bin/env nextflow

/*
 * VCF to VAF Bedgraph conversion module
 */

nextflow.enable.dsl=2

process VCF2BEDGRAPH_VAF {
  /*
  Convert VCF file to VAF (Variant Allele Frequency) bedgraph format

  Parameters
  ----------
  barcode : str
    Sample barcode/ID
  vcf : path
    Input VCF file from DeepVariant
  vcf_index : path
    VCF index file (.tbi)

  Returns
  -------
  Tuple of barcode, VAF bedgraph file, and its index
  */

  tag "$barcode"

  publishDir "${params.data}/samples/${barcode}/sequences",
    mode: 'copy',
    pattern: "${barcode}.vaf.bedgraph.gz*"

  container 'docker://ghcr.io/nvnieuwk/vcf2bedgraph:latest'
  
  label 'vcf2bedgraph'

  input:
  tuple val(barcode), path(vcf), path(vcf_index)

  output:
  tuple val(barcode), path("${output_file}.gz"), path("${output_file}.gz.tbi"), emit: vaf_bedgraph

  script:
  output_file = "${barcode}.vaf.bedgraph"

  """
  # Convert VCF to VAF bedgraph
  vcf2bedgraph \\
      ${vcf} \\
      --output ${output_file} 
  """

  stub:
  output_file = "${barcode}.vaf.bedgraph"
  """
  touch ${output_file}.gz
  touch ${output_file}.gz.tbi
  """
}
