#!/usr/bin/env nextflow

/*
 * WisecondorX cohort merge module
 * Merges all family WisecondorX aberrations BED files into cohort-level BED file
 */

nextflow.enable.dsl=2

process MERGE_COHORT_ABERRATIONS {
  /*
  Merge all family WisecondorX aberrations BED files into cohort-level BED file

  Parameters
  ----------
  cohort_name : val
    Name of the cohort
  family_bed_files : list
    List of family aberrations BED files

  Returns
  -------
  Tuple of cohort name and merged BED file
  */

  tag "${cohort_name}"

  publishDir "${params.data}/cohorts/${cohort_name}/svs/wisecondorx",
    mode: 'copy',
    pattern: "${cohort_name}_aberrations.bed"

  container 'docker://ubuntu:22.04'
  
  label 'merge_cohort_aberrations'

  input:
  tuple val(cohort_name), path(family_bed_files)

  output:
  tuple val(cohort_name), path("${output_bed}"), emit: cohort_aberrations_bed

  script:
  output_bed = "${cohort_name}_aberrations.bed"

  """
  # List all BED files (sorted for consistency)
  ls -1 *_aberrations.annotated.bed | sort -u > file_list.txt
  
  # Get the first file to check for header
  first_file=\$(head -n 1 file_list.txt)
  first_line=\$(head -n 1 "\${first_file}")
  
  # Extract first field (column) from the first line
  first_field=\$(echo "\${first_line}" | awk '{print \$1}')
  
  # Check if first field is a typical header indicator (case-insensitive)
  # Headers typically start with: #, track, browser, chr, chrom, chromosome (not a number or X/Y)
  if echo "\${first_field}" | grep -qiE '^(#|track|browser|chr|chrom|chromosome)\$'; then
    # Has header - write it with added family_id column
    echo -e "\${first_line}\tfamily_id" > ${output_bed}
    
    # Concatenate all files, skipping their first line (header) and adding family_id column
    while IFS= read -r file; do
      # Extract family ID from filename (remove _aberrations.annotated.bed suffix)
      fid=\$(basename "\${file}" _aberrations.annotated.bed)
      # Add family_id column to each line
      tail -n +2 "\${file}" | awk -v fid="\${fid}" '{print \$0 "\t" fid}' >> ${output_bed}
    done < file_list.txt
    
    # Sort only the data lines (skip header at line 1)
    head -n 1 ${output_bed} > ${output_bed}.tmp
    tail -n +2 ${output_bed} | sort -k1,1V -k2,2n >> ${output_bed}.tmp
    mv ${output_bed}.tmp ${output_bed}
  else
    # No header detected, just concatenate with family_id column and sort
    > ${output_bed}
    while IFS= read -r file; do
      # Extract family ID from filename
      fid=\$(basename "\${file}" _aberrations.annotated.bed)
      # Add family_id column to each line
      awk -v fid="\${fid}" '{print \$0 "\t" fid}' "\${file}" >> ${output_bed}
    done < file_list.txt
    
    # Sort the file
    sort -k1,1V -k2,2n ${output_bed} > ${output_bed}.tmp
    mv ${output_bed}.tmp ${output_bed}
  fi
  """

  stub:
  output_bed = "${cohort_name}_aberrations.bed"
  """
  echo "chr\tstart\tend" > ${output_bed}
  """
}
