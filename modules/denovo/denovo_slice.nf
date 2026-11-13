#!/usr/bin/env nextflow

/*
 * De novo: per-family read slices for parent of origin
 * Each child's reads around its candidates, from its CRAM to an indexed BAM
 */

nextflow.enable.dsl=2

process DENOVO_SLICE {
  /*
  Copy each child's reads around its candidates to a BAM, for Unfazed

  Unfazed's biocontainer has pysam 0.16 (htslib 1.10), which cannot read CRAM 3.1, and it reopens
  the alignment file for every candidate: on the CRAMs of CANDY_mpx_ghent's first run it read no
  read at all, and ran out of 4 GB on half the families. samtools copies, for every child of
  kids.tsv, the reads of its <kid>.windows.bed (6 kb around its candidates, from DENOVO_EVIDENCE)
  into slices/<kid>.slice.bam, indexed: a few MB to a few hundred MB per child.

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

  Returns
  -------
  Tuple of family ID and the slices directory
  */

  tag "${fid}"

  label 'denovo_slice'

  input:
  tuple val(fid), path(unfazed, stageAs: 'unfazed'), path(crams), val(label)

  output:
  tuple val(fid), path('slices'), emit: slices

  script:
  """
  set -euo pipefail

  mkdir -p slices
  while read -r kid; do
    [[ -n "\${kid}" ]] || continue
    crams=( "\${kid}".*.cram )
    if [[ \${#crams[@]} -ne 1 || ! -e "\${crams[0]}" ]]; then
      echo "ERROR: expected one CRAM named \${kid}.<ref_name>.cram" >&2
      exit 1
    fi
    samtools view -b -M -L "unfazed/\${kid}.windows.bed" -T ${params.ref} -@ ${task.cpus} \\
      -o "slices/\${kid}.slice.bam" "\${crams[0]}"
    samtools index "slices/\${kid}.slice.bam"
    echo "\${kid}: \$(samtools view -c "slices/\${kid}.slice.bam") reads" >&2
  done < unfazed/kids.tsv
  """

  stub:
  """
  mkdir -p slices
  """
}
