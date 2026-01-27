#!/usr/bin/env nextflow

/*
 * WisecondorX aberrations annotation module
 * Annotates family aberrations BED files with gene and exon information from gencode
 */

nextflow.enable.dsl=2

process ANNOTATE_ABERRATIONS {
  /*
  Annotate family WisecondorX aberrations BED files with gencode gene and exon information

  Parameters
  ----------
  fid : val
    Family ID
  aberrations_bed : path
    Family aberrations BED file
  annotation_path : val
    Directory containing gencode annotation files
  annotation_gencode : val
    Gencode version identifier (e.g., gencode.v47.basic)

  Returns
  -------
  Tuple of family ID and annotated BED file
  */

  tag "${fid}"

  publishDir "${params.data}/families/${fid}/svs/wisecondorx",
    mode: 'copy',
    pattern: "${fid}_aberrations.annotated.bed"

  container 'docker://staphb/bedtools:latest'
  
  label 'annotate_aberrations'

  input:
  tuple val(fid), path(aberrations_bed), val(annotation_path), val(annotation_gencode)

  output:
  tuple val(fid), path("${output_bed}"), emit: annotated_aberrations

  script:
  output_bed = "${fid}_aberrations.annotated.bed"
  gencode_dir = "${annotation_path}/gencode"
  symbol_genes = "${gencode_dir}/${annotation_gencode}.symbol.genes.bed.gz"
  symbol_exons = "${gencode_dir}/${annotation_gencode}.symbol.exons.bed.gz"
  ensg_genes = "${gencode_dir}/${annotation_gencode}.ensg.genes.bed.gz"
  ensg_exons = "${gencode_dir}/${annotation_gencode}.ensg.exons.bed.gz"

  """
  #!/bin/bash
  set -euo pipefail
  
  # Extract header if present
  first_field=\$(head -n 1 ${aberrations_bed} | awk '{print \$1}')
  
  if echo "\${first_field}" | grep -qiE '^(#|track|browser|chr|chrom|chromosome)\$'; then
    HAS_HEADER=true
    # Save header and add new columns
    head -n 1 ${aberrations_bed} | awk '{print \$0 "\\tgenic_symbol\\tgenic_ensg\\texonic_symbol\\texonic_ensg"}' > header.txt
    # Extract data (skip header)
    tail -n +2 ${aberrations_bed} > data.bed
  else
    HAS_HEADER=false
    cp ${aberrations_bed} data.bed
  fi
  
  # Check if data.bed is empty
  if [ ! -s data.bed ]; then
    # Empty file - just create output with header if needed
    if [ "\${HAS_HEADER}" = "true" ]; then
      cp header.txt ${output_bed}
    else
      touch ${output_bed}
    fi
    exit 0
  fi
  
  # Add line numbers to track original order
  awk '{print \$0 "\\t" NR}' data.bed > data_numbered.bed
  
  # Intersect with each annotation file and collect results
  # 1. Genic symbols
  bedtools intersect -a data_numbered.bed -b ${symbol_genes} -wa -wb | \
    awk '{line_num=\$(NF); gene=\$(NF-3); print line_num "\\t" gene}' | \
    sort -u | \
    sort -k1,1n | \
    awk '{genes[NR]=\$1; symbols[NR]=symbols[NR] (symbols[NR] ? "," : "") \$2} 
         END {for(i=1; i<=NR; i++) if(genes[i]!=genes[i-1] || i==1) print genes[i] "\\t" symbols[i]}' > genic_symbol.txt
  
  # 2. Genic ENSG
  bedtools intersect -a data_numbered.bed -b ${ensg_genes} -wa -wb | \
    awk '{line_num=\$(NF); gene=\$(NF-3); print line_num "\\t" gene}' | \
    sort -u | \
    sort -k1,1n | \
    awk '{genes[NR]=\$1; symbols[NR]=symbols[NR] (symbols[NR] ? "," : "") \$2} 
         END {for(i=1; i<=NR; i++) if(genes[i]!=genes[i-1] || i==1) print genes[i] "\\t" symbols[i]}' > genic_ensg.txt
  
  # 3. Exonic symbols
  bedtools intersect -a data_numbered.bed -b ${symbol_exons} -wa -wb | \
    awk '{line_num=\$(NF); gene=\$(NF-3); print line_num "\\t" gene}' | \
    sort -u | \
    sort -k1,1n | \
    awk '{genes[NR]=\$1; symbols[NR]=symbols[NR] (symbols[NR] ? "," : "") \$2} 
         END {for(i=1; i<=NR; i++) if(genes[i]!=genes[i-1] || i==1) print genes[i] "\\t" symbols[i]}' > exonic_symbol.txt
  
  # 4. Exonic ENSG
  bedtools intersect -a data_numbered.bed -b ${ensg_exons} -wa -wb | \
    awk '{line_num=\$(NF); gene=\$(NF-3); print line_num "\\t" gene}' | \
    sort -u | \
    sort -k1,1n | \
    awk '{genes[NR]=\$1; symbols[NR]=symbols[NR] (symbols[NR] ? "," : "") \$2} 
         END {for(i=1; i<=NR; i++) if(genes[i]!=genes[i-1] || i==1) print genes[i] "\\t" symbols[i]}' > exonic_ensg.txt
  
  # Merge annotations with original data
  awk '{print NR "\\t" \$0}' data.bed | \
    sort -k1,1n | \
    awk 'BEGIN {
      # Read annotations
      while((getline < "genic_symbol.txt") > 0) genic_sym[\$1]=\$2;
      while((getline < "genic_ensg.txt") > 0) genic_ensg[\$1]=\$2;
      while((getline < "exonic_symbol.txt") > 0) exonic_sym[\$1]=\$2;
      while((getline < "exonic_ensg.txt") > 0) exonic_ensg[\$1]=\$2;
    }
    {
      line_num=\$1;
      \$1="";
      sub(/^[ \\t]+/, "");
      print \$0 "\\t" (genic_sym[line_num]?genic_sym[line_num]:"") "\\t" \
                     (genic_ensg[line_num]?genic_ensg[line_num]:"") "\\t" \
                     (exonic_sym[line_num]?exonic_sym[line_num]:"") "\\t" \
                     (exonic_ensg[line_num]?exonic_ensg[line_num]:"");
    }' > annotated_data.bed
  
  # Combine header and data
  if [ "\${HAS_HEADER}" = "true" ]; then
    cat header.txt annotated_data.bed > ${output_bed}
  else
    mv annotated_data.bed ${output_bed}
  fi
  """

  stub:
  output_bed = "${fid}_aberrations.annotated.bed"
  """
  echo "chr\tstart\tend\tgenic_symbol\tgenic_ensg\texonic_symbol\texonic_ensg" > ${output_bed}
  """
}
