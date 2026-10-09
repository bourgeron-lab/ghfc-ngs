#!/usr/bin/env nextflow

/*
 * ROH: per-family runs of homozygosity and inbreeding module
 * Calls ROH on the family's ancestry panel genotypes and scores them against the reference
 */

nextflow.enable.dsl=2

process ROH_CALL {
  /*
  Call runs of homozygosity and score each family member's inbreeding

  Reads the family panel BCF the ancestry step built from gVCFs, so a covered site
  with no variant is a real 0/0. The family's common*.bcf would not do: it holds
  only sites where someone in the family carries an alt allele, which hides most of
  the homozygous evidence and makes F_ROH depend on family size. ancestry-pgs
  refuses such an input outright.

  Each member is called with the allele frequencies of its own ancestry region (from
  the family's ancestry table) and compared with unrelated reference individuals of
  that region. The cohort pedigree is staged whole: the tool only reads the rows of
  the family's own samples, and needs no family pedigree from deepvariant_family.

  Parameters
  ----------
  fid : val
    Family ID
  bcf : path
    Family panel genotype BCF file (from the ancestry step)
  bcf_index : path
    BCF index file (.csi)
  ancestry_tsv : path
    The family's ancestry table from ANCESTRY_PCS
  pedigree : path
    Cohort pedigree file
  reference : val
    Path to the mounted ancestry-pgs reference bundle, with its roh/ tier
  panel_name : val
    Label identifying this panel/threshold combination in output file names
  genes_bed : val
    Gene symbol BED, or empty to skip the gene table
  gene_ids_bed : val
    Gene ID BED with the same intervals, or empty
  min_gene_roh_mb : val
    Minimum ROH length intersected with genes

  Returns
  -------
  Tuple of family ID and the per-sample, segment, gene and kinship tables, and the QC report
  */

  tag "$fid"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/roh"
  },
    mode: 'copy',
    pattern: "${fid}.${panel_name}.{froh.tsv,roh.tsv,roh_genes.tsv,kinship.tsv,roh.qc.json}"

  label 'roh_call'

  input:
  tuple val(fid), path(bcf), path(bcf_index), path(ancestry_tsv), path(pedigree), val(reference), val(panel_name), val(genes_bed), val(gene_ids_bed), val(min_gene_roh_mb)

  output:
  tuple val(fid), path("${out_prefix}.froh.tsv"), emit: froh
  tuple val(fid), path("${out_prefix}.roh.tsv"), emit: roh
  tuple val(fid), path("${out_prefix}.roh_genes.tsv"), emit: roh_genes, optional: true
  // Written only when the family has two or more members
  tuple val(fid), path("${out_prefix}.kinship.tsv"), emit: kinship, optional: true
  tuple val(fid), path("${out_prefix}.roh.qc.json"), emit: qc

  script:
  out_prefix = "${fid}.${panel_name}"
  def genes_args = genes_bed ? "--genes ${genes_bed} --min-gene-roh-mb ${min_gene_roh_mb}" : ""
  def gene_ids_args = (genes_bed && gene_ids_bed) ? "--gene-ids ${gene_ids_bed}" : ""

  """
  set -euo pipefail

  ancestry-pgs roh \\
      --bcf ${bcf} \\
      --reference ${reference} \\
      --ancestry ${ancestry_tsv} \\
      --pedigree ${pedigree} \\
      ${genes_args} ${gene_ids_args} \\
      --out-prefix ${out_prefix} \\
      --threads ${task.cpus}
  """

  stub:
  out_prefix = "${fid}.${panel_name}"
  def genes_stub = genes_bed ? "printf 'IID\\tchrom\\tgene\\n' > ${out_prefix}.roh_genes.tsv" : "true"
  """
  printf 'IID\\tregion\\tfroh_1.5mb\\tfroh_category\\tfather\\tmother\\n${fid}_s1\\tEUR\\t0.001\\tnone\\t.\\t.\\n' > ${out_prefix}.froh.tsv
  printf 'IID\\tchrom\\tstart\\tend\\tlength_mb\\n' > ${out_prefix}.roh.tsv
  ${genes_stub}
  echo '{}' > ${out_prefix}.roh.qc.json
  """
}
