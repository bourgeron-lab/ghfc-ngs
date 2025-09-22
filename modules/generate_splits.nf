#!/usr/bin/env nextflow

/*
 * Script to generate split files for VCF chunking
 */

nextflow.enable.dsl=2

process generateSplits {
  /*
  Generate split files that each contain bin_size number of variants

  Returns
  -------
  Tuple of VCF, VCF index, split files, vep config file
  */

  cpus params.cpus
  label 'bcftools'

  input:
  tuple val(meta), path(vcf), path(vcf_index), path(vep_config)

  output:
  tuple val(meta), path(vcf), path(vcf_index), path("x*"), path(vep_config)

  shell:
  """
  bcftools query -f'%CHROM\t%POS\n' ${vcf} | uniq | split -a 3 -l ${params.bin_size}
  """
}