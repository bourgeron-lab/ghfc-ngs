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
  # List all BED files (sorted for consistency)
  ls -1 *_aberrations.bed | sort -u > file_list.txt
  
  # Get the first file to extract the header (if it exists)
  first_file=\$(head -n 1 file_list.txt)
  
  # Check if the first line is a header
  # Headers can be: comment lines (#, track, browser) or column names (chr, chrom, start, end, etc.)
  first_line=\$(head -n 1 "\${first_file}")
  
  # Check if it's a comment-style header OR starts with typical BED column names
  if echo "\${first_line}" | grep -qiE '^(#|track|browser|chr\s|chrom\s|chromosome\s)'; then
    # Has header - write it once from first file
    echo "\${first_line}" > ${output_bed}
    
    # Concatenate all files, skipping their first line (header)
    while IFS= read -r file; do
      tail -n +2 "\${file}" >> ${output_bed}
    done < file_list.txt
  else
    # No header detected, just concatenate all files
    cat *_aberrations.bed > ${output_bed}
  fi
  
  # Sort the output by chromosome and position (excluding header line if present)
  if head -n 1 ${output_bed} | grep -qiE '^(#|track|browser|chr\s|chrom\s|chromosome\s)'; then
    # Extract header
    head -n 1 ${output_bed} > ${output_bed}.tmp
    # Sort data lines and append
    tail -n +2 ${output_bed} | sort -k1,1 -k2,2n >> ${output_bed}.tmp
    mv ${output_bed}.tmp ${output_bed}
  else
    # No header, just sort
    sort -k1,1 -k2,2n ${output_bed} > ${output_bed}.tmp
    mv ${output_bed}.tmp ${output_bed}
  fi
  """

  stub:
  output_bed = "${fid}_aberrations.bed"
  """
  echo "chr\tstart\tend" > ${output_bed}
  """
}
