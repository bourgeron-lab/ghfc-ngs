#!/usr/bin/env nextflow

/*
 * UPD: per-family module
 * The uniparental disomy segments of every member of a family, from its norm.bcf
 */

nextflow.enable.dsl=2

process UPD_FAMILY {
  /*
  Find the UPD segments of a family's members

  The snapshot, built by the plan from the pedigree (lib/InputsSnapshot.groovy), is written
  as <prefix>.inputs.tsv: the family's rows, normalised, which members are in the BCF and
  which have a coverage JSON, and which roh tables are read. The plan rebuilds the family when
  the one it builds for a run differs from the one on disk, so a pedigree edit reruns it.

  The family's biallelic SNVs on the autosomes and chrX get their gnomAD frequency (the
  upd_af_field of gnomad_file, as INFO/UPD_AF), and stream into upd_segments family, which drops
  those outside the roh panel's covered segments when there is a mask: an HMM per member with
  both parents in the BCF (trio) or one (duo), and the roh step's ROH segments for a member with
  none (solo). See its docstring for the states, the score and the rules. It also writes the
  windows table, the sites of each member's model counted in 500 kb windows, for genetics-viz to
  draw the evidence behind the calls; the cohort does not gather it.

  Parameters
  ----------
  fid : val
    Family ID
  norm : path
    The family's norm.bcf and its .csi
  coverage : path
    The coverage JSONs of the members with one; may be empty
  roh : path
    The family's roh.tsv and froh.tsv from the roh step, or nothing
  snapshot : val
    The family's inputs snapshot, as text
  gnomad_file : val
    The gnomAD BCF, read in place as the annotation step does
  mask : val
    The roh panel's sites (roh/sites.tsv.gz of ancestry_reference), read in place; sites
    outside its covered segments are dropped. Empty for no mask
  label : val
    upd_label, identifying the analysis in output file names
  segments : path
    upd_segments, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of family ID and its segments, events, members and windows tables, and the inputs
  snapshot
  */

  tag "${fid}"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/upd"
  },
    mode: 'copy',
    pattern: "${fid}.${label}.{segments,events,members,inputs,windows}.tsv"

  label 'upd_family'

  input:
  // Staged into subdirectories, so a family named after one of its members cannot have an
  // input and an output under the same name
  tuple val(fid), path(norm, stageAs: 'bcf/*'), path(coverage, stageAs: 'coverage/*'), path(roh, stageAs: 'roh/*'),
        val(snapshot), val(gnomad_file), val(mask), val(label), path(segments)

  output:
  tuple val(fid), path("${prefix}.segments.tsv"), path("${prefix}.events.tsv"), path("${prefix}.members.tsv"),
        path("${prefix}.windows.tsv"), path("${prefix}.inputs.tsv"), emit: tables

  script:
  prefix = "${fid}.${label}"
  def bcf = (norm instanceof List ? norm : [norm]).find { staged -> staged.name.endsWith('.bcf') }
  def roh_files = roh instanceof List ? roh : (roh ? [roh] : [])
  def roh_tsv = roh_files.find { staged -> staged.name.endsWith('.roh.tsv') }
  def froh_tsv = roh_files.find { staged -> staged.name.endsWith('.froh.tsv') }
  def field = params.upd_af_field

  """
  set -euo pipefail

  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}

  # A field the gnomAD file does not have would annotate nothing, and every site would get the
  # frequency floor without a word
  if ! bcftools view -h ${gnomad_file} | grep -q '^##INFO=<ID=${field},'; then
    echo "ERROR: ${gnomad_file} has no INFO/${field} (upd_af_field)" >&2
    exit 1
  fi

  roh=()
  ${roh_tsv && froh_tsv ? "roh=(--roh ${roh_tsv} --froh ${froh_tsv})" : ''}

  # Annotated before the filter: bcftools annotate needs an indexed target
  regions=\$(bcftools index -s ${bcf} | cut -f1 | grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X)\$' | paste -sd, -)
  bcftools annotate -a ${gnomad_file} -c INFO/UPD_AF:=INFO/${field} -r "\${regions}" -Ou ${bcf} \\
    | bcftools view -v snps -m2 -M2 -Ou \\
    | bcftools query -f '%CHROM\\t%POS\\t%INFO/UPD_AF[\\t%GT\\t%GQ\\t%DP]\\n' \\
    | python3 ${segments} family \\
        --snapshot ${prefix}.inputs.tsv \\
        --samples "\$(bcftools query -l ${bcf} | paste -sd, -)" \\
        --coverage-dir coverage \\
        ${mask ? "--mask ${mask}" : ''} \\
        \${roh[@]+"\${roh[@]}"} \\
        --out-prefix ${prefix}
  """

  stub:
  prefix = "${fid}.${label}"
  """
  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}
  printf 'family_id\\tsample\\tchrom\\tstart\\tend\\n' > ${prefix}.segments.tsv
  printf 'family_id\\tsample\\tchrom\\tstart\\tend\\n' > ${prefix}.events.tsv
  printf 'family_id\\tsample\\tmode\\tstatus\\n${fid}\\t${fid}_s1\\ttrio\\tanalysed\\n' > ${prefix}.members.tsv
  printf 'family_id\\tsample\\tmode\\tchrom\\tstart\\tend\\tn_sites\\tn_het\\tn_hom_alt\\tn_mat\\tn_pat\\tn_lacks_father\\tn_lacks_mother\\tdepth_ratio\\n' > ${prefix}.windows.tsv
  """
}
