#!/usr/bin/env nextflow

/*
 * WisecondorX aberrations chromosome reformat module
 * Adds "chr" prefix to chromosome names in aberrations BED files
 */

nextflow.enable.dsl=2

process REFORMAT_CHR {
  /*
  Reformat aberrations BED files to add "chr" prefix to chromosome names

  Parameters
  ----------
  barcode : val
    Sample barcode
  aberrations_bed : path
    Individual aberrations BED file

  Returns
  -------
  Tuple of barcode and reformatted BED file
  */

  tag "${barcode}"

  publishDir "${params.data}/samples/${barcode}/svs/wisecondorx",
    mode: 'copy',
    pattern: "${barcode}_aberrations.chr.bed"

  container 'docker://ubuntu:22.04'
  
  label 'reformat_chr'

  input:
  tuple val(barcode), path(aberrations_bed)

  output:
  tuple val(barcode), path("${output_bed}"), emit: chr_aberrations

  script:
  output_bed = "${barcode}_aberrations.chr.bed"

  """
  #!/bin/bash
  set -euo pipefail
  
  # Extract header if present
  first_field=\$(head -n 1 ${aberrations_bed} | awk '{print \$1}')
  
  if echo "\${first_field}" | grep -qiE '^(#|track|browser|chr|chrom|chromosome)\$'; then
    # Has header - preserve it as-is
    head -n 1 ${aberrations_bed} > ${output_bed}
    
    # Add chr prefix to data lines (skip header)
    tail -n +2 ${aberrations_bed} | awk 'BEGIN{OFS="\\t"} {
      # Only add chr prefix if not already present
      if (\$1 !~ /^chr/) {
        \$1 = "chr" \$1
      }
      print
    }' >> ${output_bed}
  else
    # No header - just add chr prefix to all lines
    awk 'BEGIN{OFS="\\t"} {
      # Only add chr prefix if not already present
      if (\$1 !~ /^chr/) {
        \$1 = "chr" \$1
      }
      print
    }' ${aberrations_bed} > ${output_bed}
  fi
  """

  stub:
  output_bed = "${barcode}_aberrations.chr.bed"
  """
  echo "chr1\t1000\t2000" > ${output_bed}
  """
}
