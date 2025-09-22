#!/usr/bin/env nextflow

/*
 * Script to split a VCF file into multiple smaller VCFs
 */

nextflow.enable.dsl=2

// defaults
prefix = "out"

process splitVCF {
  /*
  Split VCF file into multiple smaller VCFs based on split files

  Returns
  -------
  Tuple of original VCF, split VCF files, split VCF index files, vep config file
  */

  label 'bcftools'

  input:
  tuple val(meta), path(vcf), path(vcf_index), path(split_file), val(vep_config)

  output:
  tuple val(meta), val("${vcf}"), path("${prefix}*.vcf.gz"), path("${prefix}*.vcf.gz.{tbi,csi}"), val(vep_config)

  afterScript 'rm x*'

  script:
  index_type = meta.index_type
  index_flag = index_type == "tbi" ? "-t" : "-c"

  """
  bcftools view --no-version -T ${split_file} -Oz ${vcf} > ${prefix}.${split_file}.vcf.gz
  bcftools index ${index_flag} ${prefix}.${split_file}.vcf.gz
  """
}