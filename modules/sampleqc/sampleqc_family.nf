#!/usr/bin/env nextflow

/*
 * Sample QC: per-family module
 * The family's relatedness against its pedigree, every member's sex and karyotype, and its
 * parent-child Mendelian counts
 */

nextflow.enable.dsl=2

process SAMPLEQC_FAMILY {
  /*
  Check a family's pedigree and its members' sex against their genotypes and coverage

  The snapshot, built by the plan from the pedigree (lib/InputsSnapshot.groovy), is written
  as <prefix>.inputs.tsv: the family's rows, normalised, which members have an extract and
  a coverage JSON, and whether the family BCF is read. The plan rebuilds the family when the
  one it builds for a run differs from the one on disk, so a pedigree edit reruns it. It is
  also the pedigree every tool here reads.

  somalier relate gives every pair of members its IBS0, IBS2 and shared heterozygotes, and
  every member its depth and X/Y counts. With the family's norm.bcf, sampleqc_classify
  mendelian counts the Mendelian errors of every trio and the opposite homozygotes of every
  parent-child pair, at the biallelic autosomal SNVs, genotypes of GQ >= 20 and DP >= 10.
  sampleqc_classify family turns them into the three tables; see its docstring for the
  classes, the flags and their thresholds.

  Parameters
  ----------
  fid : val
    Family ID
  extracts : path
    The .somalier files of the members with one
  coverage : path
    The coverage JSONs of the members with one; may be empty
  norm : path
    The family's norm.bcf and its .csi, or nothing when the snapshot does not use it
  snapshot : val
    The family's inputs snapshot, as text
  sites : path
    somalier's sites VCF
  label : val
    sampleqc_label, identifying the QC in output file names
  classify : path
    sampleqc_classify, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of family ID and its samples, pairs and flags tables, and the inputs snapshot
  */

  tag "${fid}"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/qc"
  },
    mode: 'copy',
    pattern: "${fid}.${label}.{samples,pairs,flags,inputs}.tsv"

  label 'sampleqc_family'

  input:
  // Staged into subdirectories, so a family named after one of its members cannot have an
  // input and an output under the same name
  tuple val(fid), path(extracts, stageAs: 'extracts/*'), path(coverage, stageAs: 'coverage/*'),
        path(norm, stageAs: 'bcf/*'), val(snapshot), path(sites), val(label), path(classify)

  output:
  tuple val(fid), path("${prefix}.samples.tsv"), path("${prefix}.pairs.tsv"), path("${prefix}.flags.tsv"),
        path("${prefix}.inputs.tsv"), emit: tables

  script:
  prefix = "${fid}.${label}"
  def bcf = (norm instanceof List ? norm : [norm]).find { staged -> staged.name.endsWith('.bcf') }

  """
  set -euo pipefail

  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}

  # somalier's pedigree: the snapshot's members with an extract
  awk -F'\\t' -v OFS='\\t' '!/^#/ && \$6 == "yes" { print \$1, \$2, \$3, \$4, \$5, -9 }' ${prefix}.inputs.tsv > family.ped
  find extracts -name '*.somalier' | sort > extracts.list
  if [ "\$(wc -l < extracts.list)" -ne "\$(wc -l < family.ped)" ]; then
    echo "ERROR: ${fid} has \$(wc -l < extracts.list | tr -d ' ') extracts staged, and \$(wc -l < family.ped | tr -d ' ') members with one in its snapshot" >&2
    exit 1
  fi
  somalier relate --sites ${sites} --ped family.ped -o relate extracts.list

  mendelian=()
  ${bcf ? """
  autosomes=\$(bcftools index -s ${bcf} | cut -f1 | grep -E '^(chr)?([1-9]|1[0-9]|2[0-2])\$' | paste -sd, -)
  bcftools view -v snps -m2 -M2 -r "\${autosomes}" -Ou ${bcf} \\
    | bcftools query -f '[%GT\\t%GQ\\t%DP\\t]\\n' \\
    | python3 ${classify} mendelian --snapshot ${prefix}.inputs.tsv \\
        --samples "\$(bcftools query -l ${bcf} | paste -sd, -)" --out mendelian.tsv
  mendelian=(--mendelian mendelian.tsv)""" : ''}

  python3 ${classify} family \\
      --snapshot ${prefix}.inputs.tsv \\
      --somalier-samples relate.samples.tsv \\
      --somalier-pairs relate.pairs.tsv \\
      --coverage-dir coverage \\
      \${mendelian[@]+"\${mendelian[@]}"} \\
      --out-prefix ${prefix}
  """

  stub:
  prefix = "${fid}.${label}"
  """
  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}
  printf 'family_id\\tsample\\tflags\\n${fid}\\t${fid}_s1\\t.\\n' > ${prefix}.samples.tsv
  printf 'family_id_a\\tsample_a\\tfamily_id_b\\tsample_b\\tstatus\\n' > ${prefix}.pairs.tsv
  printf 'level\\tfamily_id\\tsample_a\\tsample_b\\tflag\\tseverity\\tevidence\\n' > ${prefix}.flags.tsv
  """
}
