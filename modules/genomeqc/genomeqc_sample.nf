#!/usr/bin/env nextflow

/*
 * Genome QC: per-sample module
 * What kind of DNA and library the sample is, from its CRAM
 */

nextflow.enable.dsl=2

process GENOMEQC_SAMPLE {
  /*
  Count a CRAM's reads on chrEBV, chrM and chr20, and take samtools stats of its reads on chr20

  The counts are of primary, non-duplicate reads of MAPQ >= 20, each an index jump to one
  contig. chr20 stands for the genome: an EBV-transformed cell line carries tens of EBV genomes
  per cell, blood almost none, and the copies per cell are 2 x (EBV reads per base) / (chr20
  reads per base). A contig the CRAM's header lacks is counted NA. samtools stats on chr20,
  about 2% of the genome, gives the read length, the duplicate, MAPQ0 and properly paired
  shares, the mismatch rate, the insert size and the indels in the reads' alignments, at a
  fraction of a whole-CRAM pass. genomeqc_tables sample turns both into one row.

  Parameters
  ----------
  barcode : val
    Sample ID
  cram : path
    Sample CRAM file
  crai : path
    CRAM index file
  ref : path
    The reference the CRAM was aligned to
  ref_fai : path
    Its index, staged beside it
  label : val
    genomeqc_label, identifying the QC in output file names
  tables_script : path
    genomeqc_tables, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of barcode, its metrics table and its samtools stats
  */

  tag "${barcode}"

  publishDir {
    def (s1, s2) = Sharding.getShards(barcode)
    "${params.data}/samples/${s1}/${s2}/${barcode}/qc/genomeqc"
  },
    mode: 'copy',
    pattern: "${barcode}.${label}.{tsv,stats.txt}"

  label 'genomeqc_sample'

  input:
  tuple val(barcode), path(cram), path(crai), path(ref), path(ref_fai), val(label), path(tables_script)

  output:
  tuple val(barcode), path("${barcode}.${label}.tsv"), path("${barcode}.${label}.stats.txt"), emit: metrics

  script:
  def prefix = "${barcode}.${label}"
  """
  set -euo pipefail

  samtools view -H --reference ${ref} ${cram} > header.sam
  : > counts.tsv
  for contig in chrEBV chrM chr20; do
    if awk -F'\\t' -v c="SN:\${contig}" '\$1 == "@SQ" && \$2 == c { found = 1 } END { exit !found }' header.sam; then
      n=\$(samtools view -c -F 0xF04 -q 20 --reference ${ref} ${cram} "\${contig}")
    else
      n=NA
    fi
    printf '%s\\t%s\\n' "\${contig}" "\${n}" >> counts.tsv
  done

  samtools stats -@ ${task.cpus} --reference ${ref} ${cram} chr20 > ${prefix}.stats.txt

  python3 ${tables_script} sample --barcode ${barcode} --counts counts.tsv --stats ${prefix}.stats.txt \\
      --out ${prefix}.tsv
  """

  stub:
  def prefix = "${barcode}.${label}"
  """
  printf 'sample\\tread_length\\tebv_copies_per_cell\\n${barcode}\\t151\\t0.001\\n' > ${prefix}.tsv
  echo stub > ${prefix}.stats.txt
  """
}
