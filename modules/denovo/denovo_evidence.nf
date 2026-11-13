#!/usr/bin/env nextflow

/*
 * De novo: per-family evidence
 * Masks, clusters and read-level evidence for every candidate of the family's Wombat de novo
 * table, the callable genome of each child, and Unfazed's inputs
 */

nextflow.enable.dsl=2

process DENOVO_EVIDENCE {
  /*
  Gather the evidence on a family's de novo candidates

  denovo_tables evidence (see its docstring for the flags and their thresholds):
    - one flag per mask of the bundle a candidate overlaps, CLUSTER for two candidates of a child
      within 20 bp, and the child's other calls within 20 bp in the norm.bcf;
    - a pileup of every member with a CRAM at the candidates, twice (reads of mapping and base
      quality 20, and every read): parental alternate reads, siblings' reads, the child's VAF,
      strand balance and low-quality support;
    - each child's callable genome from the coverage bedgraphs, less the masks' union;
    - Unfazed's inputs, for DENOVO_ORIGIN.
  Nothing is published: DENOVO_FAMILY writes the family's tables.

  Parameters
  ----------
  fid : val
    Family ID
  table : path
    The family's Wombat de novo table (denovo_wombat_config)
  norm : path
    The family's norm.bcf and its .csi
  crams : path
    The CRAMs and .crai of the members with one, under their own names
  bedgraphs : path
    The coverage bedgraphs of the members with one
  snapshot : val
    The family's inputs snapshot, as text: its pedigree rows, with which members have a CRAM
    and a bedgraph (lib/InputsSnapshot.groovy)
  masks : path
    The mask bundle directory (denovo_masks)
  label : val
    denovo_label, identifying the analysis in output file names
  tables_script : path
    denovo_tables, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of family ID, its evidence and callable tables, Unfazed's input directory and the
  inputs snapshot
  */

  tag "${fid}"

  label 'denovo_evidence'

  input:
  // Staged into subdirectories, so a family named after one of its members cannot have an
  // input and an output under the same name. The CRAMs and bedgraphs stay in the working
  // directory, where denovo_tables finds them by their members' names.
  tuple val(fid), path(table, stageAs: 'wombat/*'), path(norm, stageAs: 'bcf/*'), path(crams), path(bedgraphs),
        val(snapshot), path(masks, stageAs: 'masks'), val(label), path(tables_script)

  output:
  tuple val(fid), path("${prefix}.evidence.tsv"), path("${prefix}.callable.tsv"), path('unfazed'),
        path("${prefix}.inputs.tsv"), emit: evidence

  script:
  prefix = "${fid}.${label}"
  def bcf = (norm instanceof List ? norm : [norm]).find { staged -> staged.name.endsWith('.bcf') }

  """
  set -euo pipefail

  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}

  python3 ${tables_script} evidence \\
    --table ${table} \\
    --bcf ${bcf} \\
    --snapshot ${prefix}.inputs.tsv \\
    --masks masks \\
    --ref ${params.ref} \\
    --out-prefix ${prefix} \\
    --unfazed-dir unfazed \\
    --workdir tmp
  rm -rf tmp
  """

  stub:
  prefix = "${fid}.${label}"
  """
  ${InputsSnapshot.decodeCommand(snapshot, "${prefix}.inputs.tsv")}
  mkdir -p unfazed
  touch unfazed/kids.tsv unfazed/dnms.bed
  printf 'family_id\\tchild\\tflags\\n' > ${prefix}.evidence.tsv
  printf 'family_id\\tchild\\tfather\\tmother\\tcallable_bp\\tcallable_unmasked_bp\\tstatus\\n' > ${prefix}.callable.tsv
  """
}
