#!/usr/bin/env nextflow

/*
 * Genome QC: per-family module
 * Each member's SNV and indel profile, from the family's norm.bcf
 */

nextflow.enable.dsl=2

process GENOMEQC_FAMILY {
  /*
  Profile each member's autosomal SNVs and indels in the family's norm.bcf

  genomeqc_tables family counts each member's non-reference genotypes of GQ >= 20 and DP >= 10:
  heterozygous and homozygous SNVs, Ti/Tv, the SNVs private to the member in the family, the
  share of heterozygous SNVs under 0.3 of the reads; the indels by size and type, and, with
  genomeqc_tr_bed, the share in tandem repeats and homopolymers. bcftools stats gives indel
  sizes for the whole family only, hence a pass of our own.

  The snapshot, built by the plan from the pedigree (lib/InputsSnapshot.groovy), is written as
  <prefix>.inputs.tsv, so that the plan can tell a pedigree edit or a new TR BED.

  Parameters
  ----------
  fid : val
    Family ID
  norm : path
    The family's norm.bcf and its .csi
  tr_bed : path
    genomeqc_tr_bed, or nothing
  snapshot : val
    The family's inputs snapshot, as text
  label : val
    genomeqc_label, identifying the QC in output file names
  tables_script : path
    genomeqc_tables

  Returns
  -------
  Tuple of family ID, its variants table and the snapshot
  */

  tag "${fid}"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/qc"
  },
    mode: 'copy',
    pattern: "${fid}.${label}.{variants,inputs}.tsv"

  label 'genomeqc_family'

  input:
  // Staged into subdirectories, so a family named after one of its members cannot have an
  // input and an output under the same name
  tuple val(fid), path(norm, stageAs: 'bcf/*'), path(tr_bed, stageAs: 'tr/*'), val(snapshot), val(label),
        path(tables_script)

  output:
  tuple val(fid), path("${prefix}.variants.tsv"), path("${prefix}.inputs.tsv"), emit: tables

  script:
  prefix = "${fid}.${label}"
  def bcf = (norm instanceof List ? norm : [norm]).find { staged -> staged.name.endsWith('.bcf') }
  def bed = (tr_bed instanceof List ? tr_bed : (tr_bed ? [tr_bed] : [])).find { staged -> !staged.name.endsWith('.tbi') && !staged.name.endsWith('.csi') }
  """
  set -euo pipefail

  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}

  python3 ${tables_script} family --family ${fid} --bcf ${bcf} ${bed ? "--tr-bed ${bed}" : ''} \\
      --out ${prefix}.variants.tsv
  """

  stub:
  prefix = "${fid}.${label}"
  """
  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}
  printf 'family_id\\tsample\\tn_snv\\n${fid}\\t${fid}_s1\\t1\\n' > ${prefix}.variants.tsv
  """
}
