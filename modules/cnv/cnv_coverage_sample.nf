#!/usr/bin/env nextflow

/*
 * CNV coverage: per-sample module
 * The sample's normalised coverage on a fixed bin grid, as a bigWig, with its summary JSON
 */

nextflow.enable.dsl=2

process CNV_COVERAGE_SAMPLE {
  /*
  Bin one sample's depth with mosdepth and write it normalised, for genetics-viz

  mosdepth -n --fast-mode --by <bin> -Q <mapq>, with its default flag filter, so that
  duplicate, secondary and QC-fail reads are left out. cnv_coverage sample then keeps
  chr1-22, chrX and chrY, divides every bin by the autosomal median depth, and writes every
  bin, zeros included, as 2 x depth / median rounded to round_to; see the script's
  docstring for the values and the JSON. The tag cov_b<bin>q<mapq> in both file names is
  what genetics-viz matches a sample to a panel by.

  This is separate from the MOSDEPTH process of the alignment, whose by<bin> bedgraph is
  left as it is.

  Parameters
  ----------
  barcode : val
    Sample ID
  cram : path
    Sample CRAM file
  crai : path
    CRAM index file
  ref : path
    The reference the CRAM was aligned to; its basename is recorded in the JSON
  ref_fai : path
    Its index, which gives the contig lengths of the grid
  bin : val
    cnv_cov_bin, the bin size in bp
  mapq : val
    cnv_cov_mapq, the minimum mapping quality
  round_to : val
    cnv_cov_round, the precision the normalised values are rounded to
  script : path
    cnv_coverage, staged in so it is reachable inside the container

  Returns
  -------
  Tuple of barcode, the bigWig and the JSON
  */

  tag "${barcode}"

  publishDir {
    def (s1, s2) = Sharding.getShards(barcode)
    "${params.data}/samples/${s1}/${s2}/${barcode}/coverage"
  },
    mode: 'copy',
    pattern: "${barcode}.cov_b${bin}q${mapq}.{bw,json}"

  label 'cnv_coverage_sample'

  input:
  tuple val(barcode), path(cram), path(crai), path(ref), path(ref_fai), val(bin), val(mapq), val(round_to), path(script)

  output:
  tuple val(barcode), path("${prefix}.bw"), path("${prefix}.json"), emit: coverage

  script:
  prefix = "${barcode}.cov_b${bin}q${mapq}"

  """
  set -euo pipefail

  echo "${barcode}: mosdepth started at \$(date)"
  mosdepth -t ${task.cpus} -n --fast-mode --by ${bin} -Q ${mapq} -f ${ref} mosdepth_out ${cram}
  echo "${barcode}: mosdepth done at \$(date)"

  python3 ${script} sample \\
      --regions mosdepth_out.regions.bed.gz \\
      --fasta ${ref} \\
      --fai ${ref_fai} \\
      --sample ${barcode} \\
      --bin ${bin} \\
      --mapq ${mapq} \\
      --round ${round_to} \\
      --mosdepth-version "\$(mosdepth --version)"
  echo "${barcode}: bigWig written at \$(date), \$(du -k ${prefix}.bw | cut -f1) kB"
  """

  stub:
  prefix = "${barcode}.cov_b${bin}q${mapq}"
  """
  echo stub > ${prefix}.bw
  printf '{"sample": "${barcode}", "tag": "cov_b${bin}q${mapq}", "chrX_ratio": 1.0, "chrY_ratio": 0.0}\\n' > ${prefix}.json
  """
}
