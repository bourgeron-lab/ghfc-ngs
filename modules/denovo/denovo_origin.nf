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
  each child's reads within 6 kb of its candidates as a BAM, from DENOVO_SLICE. Never the CRAMs:
  its pysam 0.16 cannot read CRAM 3.1, and on the cohort's CRAMs it read no read at all (only the
  sons' chrX and chrY candidates, which it phases by sex, came out) and ran out of 4 GB on half
  the families. A son's chrX and chrY candidates are phased by his sex (-g 38). About a quarter of
  the candidates of a 30x genome are phased; the others are unphased.

  A failure of Unfazed is written in the table (origin error for the family's children), never a
  failed task: parent of origin is a bonus.

  Parameters
  ----------
  fid : val
    Family ID
  unfazed : path
    Unfazed's input directory, from DENOVO_EVIDENCE
  slices : path
    The children's BAM slices, from DENOVO_SLICE
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
  tuple val(fid), path(unfazed, stageAs: 'unfazed'), path(slices, stageAs: 'slices'), val(label), path(tables_script)

  output:
  tuple val(fid), path("${prefix}.origin.tsv"), emit: origin

  script:
  prefix = "${fid}.${label}"
  """
  set -euo pipefail

  python3 ${tables_script} origin \\
    --family ${fid} \\
    --unfazed-dir unfazed \\
    --slices-dir slices \\
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
