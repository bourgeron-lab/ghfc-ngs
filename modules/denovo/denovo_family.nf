#!/usr/bin/env nextflow

/*
 * De novo: per-family tables
 * The family's candidates with their evidence, flags and origin, and each child's QC row
 */

nextflow.enable.dsl=2

process DENOVO_FAMILY {
  /*
  Write a family's de novo tables

  denovo_tables family joins the evidence, the origin and the callable genome:
    <prefix>.dnm.tsv  one row per child and candidate, with its flags, strict when it has none
                      but SIB_SHARED, and its parent of origin
    <prefix>.qc.tsv   one row per child with both parents sequenced: candidates, strict ones
                      (SNVs, indels, transition fraction, median VAF), parental and mosaic flags,
                      origins, and the callable unmasked genome with the rate per Gb
  The inputs snapshot is published beside them: the plan rebuilds the family when the one it
  builds for a run differs, so a pedigree edit or a new CRAM reruns it.

  Parameters
  ----------
  fid : val
    Family ID
  evidence, callable, unfazed, inputs : path
    From DENOVO_EVIDENCE
  origin : path
    From DENOVO_ORIGIN
  label : val
    denovo_label
  tables_script : path
    denovo_tables

  Returns
  -------
  Tuple of family ID, its dnm and qc tables, and the inputs snapshot
  */

  tag "${fid}"

  publishDir {
    def (s1, s2) = Sharding.getShards(fid)
    "${params.data}/families/${s1}/${s2}/${fid}/denovo"
  },
    mode: 'copy',
    pattern: "${fid}.${label}.{dnm,qc,inputs}.tsv"

  label 'denovo_family'

  input:
  tuple val(fid), path(evidence, stageAs: 'evidence/*'), path(callable, stageAs: 'callable/*'),
        path(unfazed, stageAs: 'unfazed'), path(inputs, stageAs: 'snapshot/*'), path(origin, stageAs: 'origin/*'),
        val(label), path(tables_script)

  output:
  tuple val(fid), path("${prefix}.dnm.tsv"), path("${prefix}.qc.tsv"), path("${prefix}.inputs.tsv"), emit: tables

  script:
  prefix = "${fid}.${label}"
  """
  set -euo pipefail

  cp ${inputs} ${prefix}.inputs.tsv
  python3 ${tables_script} family \\
    --evidence ${evidence} \\
    --origin ${origin} \\
    --callable ${callable} \\
    --snapshot ${prefix}.inputs.tsv \\
    --unfazed-dir unfazed \\
    --out-prefix ${prefix}
  """

  stub:
  prefix = "${fid}.${label}"
  """
  cp ${inputs} ${prefix}.inputs.tsv
  printf 'family_id\\tchild\\tchrom\\tpos\\tref\\talt\\tflags\\tstrict\\n' > ${prefix}.dnm.tsv
  printf 'family_id\\tchild\\tn_candidates\\tn_strict\\n' > ${prefix}.qc.tsv
  """
}
