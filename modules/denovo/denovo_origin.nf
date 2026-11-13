#!/usr/bin/env nextflow

/*
 * De novo: per-family parent of origin
 * Unfazed's read-backed phasing of the candidates of the children with a CRAM and both parents
 * in the family BCF
 */

nextflow.enable.dsl=2

process DENOVO_ORIGIN {
  /*
  Phase a family's de novo candidates to the parent they arose in

  Unfazed 1.0.2 (Belyeu 2021) on the files DENOVO_EVIDENCE wrote: the candidates as a BED, the
  family's biallelic SNVs within 5 kb of them as the informative sites, a PED of the trios, and
  each child's CRAM. A son's chrX and chrY candidates are phased by his sex (-g 38). About a
  quarter to a third of the candidates of a 30x genome are phased; the others are unphased.

  A failure of Unfazed is written in the table (origin error for the family's children), never a
  failed task: parent of origin is a bonus. Unfazed's pysam 0.16 reads CRAM 3.0, not 3.1.

  Parameters
  ----------
  fid : val
    Family ID
  unfazed : path
    Unfazed's input directory, from DENOVO_EVIDENCE
  crams : path
    The CRAMs and .crai of the family's members with one
  label : val
    denovo_label
  tables_script : path
    denovo_tables, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of family ID and its origin table
  */

  tag "${fid}"

  label 'denovo_origin'

  input:
  tuple val(fid), path(unfazed, stageAs: 'unfazed'), path(crams), val(label), path(tables_script)

  output:
  tuple val(fid), path("${prefix}.origin.tsv"), emit: origin

  script:
  prefix = "${fid}.${label}"
  """
  set -euo pipefail

  python3 ${tables_script} origin \\
    --family ${fid} \\
    --unfazed-dir unfazed \\
    --ref ${params.ref} \\
    --threads ${task.cpus} \\
    --out ${prefix}.origin.tsv
  """

  stub:
  prefix = "${fid}.${label}"
  """
  printf 'family_id\\tchild\\tchrom\\tstart\\tend\\torigin\\n' > ${prefix}.origin.tsv
  """
}
