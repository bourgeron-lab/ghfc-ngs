#!/usr/bin/env nextflow

/*
 * WisecondorX family merge module
 * Merges all individual WisecondorX aberrations BED files for a family
 */

nextflow.enable.dsl=2

process MERGE_FAMILY_ABERRATIONS {
  /*
  Merge all individual WisecondorX aberrations BED files into family-level BED files

  Parameters
  ----------
  fid : val
    Family ID
  barcode_list : val
    List of barcodes in the family
  bed_files : list
    List of individual aberrations BED files

  Returns
  -------
  Tuple of family ID and merged BED file
  */

  tag "${fid}"

  publishDir "${params.data}/families/${fid}/svs/wisecondorx",
    mode: 'copy',
    pattern: "${fid}_aberrations.bed"

  container 'docker://ubuntu:22.04'
  
  label 'merge_family_aberrations'

  input:
  tuple val(fid), val(barcode_list), path(bed_files)

  output:
  tuple val(fid), path("${output_bed}"), emit: family_aberrations_bed

  script:
  output_bed = "${fid}_aberrations.bed"

  """
  # Set temp directory to current work directory and memory buffer for sort
  export TMPDIR=\$(pwd)
  SORT_OPTS="-T \${TMPDIR} -S 2G"

  # List all BED files (sorted for consistency)
  ls -1 *_aberrations.chr.bed | sort -u > file_list.txt

  # Get the first file to check for header
  first_file=\$(head -n 1 file_list.txt)
  first_line=\$(head -n 1 "\${first_file}")

  # Extract first field (column) from the first line
  first_field=\$(echo "\${first_line}" | awk '{print \$1}')

  # Check if first field is a typical header indicator (case-insensitive)
  # Headers typically start with: #, track, browser, chr, chrom, chromosome (not a number or X/Y)
  if echo "\${first_field}" | grep -qiE '^(#|track|browser|chr|chrom|chromosome)\$'; then
    # Has header - write it with added barcode column
    echo -e "\${first_line}\tbarcode" > ${output_bed}

    # Concatenate all files, skipping their first line (header) and adding barcode column
    while IFS= read -r file; do
      # Extract barcode from filename (remove _aberrations.chr.bed suffix)
      barcode=\$(basename "\${file}" _aberrations.chr.bed)
      # Add barcode column to each line
      tail -n +2 "\${file}" | awk -v bc="\${barcode}" '{print \$0 "\t" bc}' >> ${output_bed}
    done < file_list.txt

    # Sort only the data lines (skip header at line 1)
    head -n 1 ${output_bed} > ${output_bed}.tmp
    tail -n +2 ${output_bed} | sort \${SORT_OPTS} -k1,1V -k2,2n >> ${output_bed}.tmp
    mv ${output_bed}.tmp ${output_bed}
  else
    # No header detected, just concatenate with barcode column and sort
    > ${output_bed}
    while IFS= read -r file; do
      # Extract barcode from filename
      barcode=\$(basename "\${file}" _aberrations.chr.bed)
      # Add barcode column to each line
      awk -v bc="\${barcode}" '{print \$0 "\t" bc}' "\${file}" >> ${output_bed}
    done < file_list.txt

    # Sort the file
    sort \${SORT_OPTS} -k1,1V -k2,2n ${output_bed} > ${output_bed}.tmp
    mv ${output_bed}.tmp ${output_bed}
  fi
  """

  stub:
  output_bed = "${fid}_aberrations.bed"
  """
  echo "chr\tstart\tend" > ${output_bed}
  """
}
