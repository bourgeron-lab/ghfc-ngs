#!/usr/bin/env nextflow

/*
 * Wombat: De novo mutation extraction module using slivar
 */

nextflow.enable.dsl=2

process DNM_EXTRACTION {
  /*
  Extract de novo mutations from family VCF files using slivar
  Filters variants based on call rate, depth, genotype quality, and variant allele frequency
  Handles sex chromosomes and pseudo-autosomal regions (PAR) appropriately

  Parameters
  ----------
  fid : str
    Family ID
  bcf : path
    Input BCF file (e.g., FID.rare.vep_config_name.annotated.bcf)
  bcf_index : path
    BCF index file (.csi)
  pedigree : path
    Family pedigree file (FID.pedigree.tsv)
  min_callrate : val
    Minimum call rate threshold (default: 0.9)
  min_dp : val
    Minimum depth (DP) threshold (default: 10)
  min_gq : val
    Minimum genotype quality (GQ) threshold (default: 20)
  min_vaf : val
    Minimum variant allele frequency (VAF) threshold (default: 0.2)
  par1_start : val
    Start position of PAR1 on chrX (default: 10001 for GRCh38)
  par1_end : val
    End position of PAR1 on chrX (default: 2781479 for GRCh38)
  par2_start : val
    Start position of PAR2 on chrX (default: 155701383 for GRCh38)
  par2_end : val
    End position of PAR2 on chrX (default: 156030895 for GRCh38)
  annotation_path : val
    Path to directory containing annotation files (e.g., LCR.bed.gz)

  Returns
  -------
  Tuple of family ID, de novo BCF file (FID.rare.vep_config_name.annotated.dnm.bcf), 
  its index, and TSV file (FID.rare.vep_config_name.annotated.dnm.tsv)
  */

  tag "$fid"

  publishDir "${params.data}/families/${fid}/vcfs",
    mode: 'copy',
    pattern: "*.dnm.{bcf,bcf.csi,tsv}"

  container 'docker://brentp/slivar:latest'
  
  label 'dnm_extraction'

  input:
  tuple val(fid), path(bcf), path(bcf_index), path(pedigree), val(min_callrate), val(min_dp), val(min_gq), val(min_vaf), val(par1_start), val(par1_end), val(par2_start), val(par2_end), val(annotation_path)

  output:
  tuple val(fid), path("${output_bcf}"), path("${output_bcf}.csi"), path("${output_tsv}"), emit: dnm_results

  script:
  // Extract base filename without .bcf extension and add .dnm
  base_name = bcf.baseName
  output_bcf = "${base_name}.dnm.bcf"
  output_tsv = "${base_name}.dnm.tsv"

  """
  # Check if pedigree file has a header (first column contains "FID")
  # If it does, remove the header and create a temp file for slivar
  if head -n 1 ${pedigree} | cut -f1 | grep -q "^FID\$"; then
    echo "Pedigree file has header, removing it for slivar"
    tail -n +2 ${pedigree} > temp_pedigree.tsv
    PED_FILE="temp_pedigree.tsv"
  else
    echo "Pedigree file has no header, using as is"
    PED_FILE="${pedigree}"
  fi

  # Extract de novo variants using slivar
  # The slivar expr command applies filters for de novo variants
  # Filters:
  # - Call rate >= ${min_callrate}
  # - DP >= ${min_dp} in child and parents
  # - GQ >= ${min_gq} in child and parents
  # - VAF (AB - allele balance) >= ${min_vaf} in child
  #
  # De novo patterns:
  # 1. Autosomal het: kid.het && mom.hom_ref && dad.hom_ref (classic de novo)
  # 2. Autosomal hom_alt (hemizygous): kid.hom_alt && mom.hom_ref && dad.hom_ref (e.g., deletion)
  # 3. X chromosome in males (non-PAR): kid.hom_alt && mom.hom_ref (kid.sex == 1)
  # 4. X chromosome in females het: kid.het && mom.hom_ref && dad.hom_ref (kid.sex == 2)
  # 5. X chromosome in females hom: kid.hom_alt && mom.hom_ref && dad.hom_ref (kid.sex == 2)
  # 6. X chromosome PAR het: kid.het && mom.hom_ref && dad.hom_ref (both sexes, PAR regions)
  # 7. X chromosome PAR hom: kid.hom_alt && mom.hom_ref && dad.hom_ref (both sexes, PAR regions)
  # 8. Y chromosome in males: kid.hom_alt && dad.hom_ref (kid.sex == 1)
  # Note: For hom_ref parents, we also enforce AB <= 0.02 to ensure true homozygous reference
  
  slivar expr \\
    --vcf ${bcf} \\
    --ped \$PED_FILE \\
    --pass-only \\
    --info "variant.call_rate >= ${min_callrate}" \\
    --out-vcf temp_output.bcf \\
    --trio "denovo_auto_het:(variant.CHROM != 'chrX' && variant.CHROM != 'chrY' && kid.het && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq} && kid.AB >= ${min_vaf} && mom.hom_ref && mom.AB <= 0.02 && mom.DP >= ${min_dp} && mom.GQ >= ${min_gq} && dad.hom_ref && dad.AB <= 0.02 && dad.DP >= ${min_dp} && dad.GQ >= ${min_gq})" \\
    --trio "denovo_auto_hom:(variant.CHROM != 'chrX' && variant.CHROM != 'chrY' && kid.hom_alt && kid.AB >= 0.98 && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq} && mom.hom_ref && mom.AB <= 0.02 && mom.DP >= ${min_dp} && mom.GQ >= ${min_gq} && dad.hom_ref && dad.AB <= 0.02 && dad.DP >= ${min_dp} && dad.GQ >= ${min_gq})" \\
    --trio "denovo_x_male:(variant.CHROM == 'chrX' && kid.sex == 1 && (variant.POS < ${par1_start} || (variant.POS > ${par1_end} && variant.POS < ${par2_start}) || variant.POS > ${par2_end}) && kid.hom_alt && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq} && mom.hom_ref && mom.AB <= 0.02 && mom.DP >= ${min_dp} && mom.GQ >= ${min_gq})" \\
    --trio "denovo_x_female_het:(variant.CHROM == 'chrX' && kid.sex == 2 && (variant.POS < ${par1_start} || (variant.POS > ${par1_end} && variant.POS < ${par2_start}) || variant.POS > ${par2_end}) && kid.het && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq} && kid.AB >= ${min_vaf} && mom.hom_ref && mom.AB <= 0.02 && mom.DP >= ${min_dp} && mom.GQ >= ${min_gq} && dad.hom_ref && dad.AB <= 0.02 && dad.DP >= ${min_dp} && dad.GQ >= ${min_gq})" \\
    --trio "denovo_x_female_hom:(variant.CHROM == 'chrX' && kid.sex == 2 && (variant.POS < ${par1_start} || (variant.POS > ${par1_end} && variant.POS < ${par2_start}) || variant.POS > ${par2_end}) && kid.hom_alt && kid.AB >= 0.98 && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq} && mom.hom_ref && mom.AB <= 0.02 && mom.DP >= ${min_dp} && mom.GQ >= ${min_gq} && dad.hom_ref && dad.AB <= 0.02 && dad.DP >= ${min_dp} && dad.GQ >= ${min_gq})" \\
    --trio "denovo_x_par_het:(variant.CHROM == 'chrX' && ((variant.POS >= ${par1_start} && variant.POS <= ${par1_end}) || (variant.POS >= ${par2_start} && variant.POS <= ${par2_end})) && kid.het && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq} && kid.AB >= ${min_vaf} && mom.hom_ref && mom.AB <= 0.02 && mom.DP >= ${min_dp} && mom.GQ >= ${min_gq} && dad.hom_ref && dad.AB <= 0.02 && dad.DP >= ${min_dp} && dad.GQ >= ${min_gq})" \\
    --trio "denovo_x_par_hom:(variant.CHROM == 'chrX' && ((variant.POS >= ${par1_start} && variant.POS <= ${par1_end}) || (variant.POS >= ${par2_start} && variant.POS <= ${par2_end})) && kid.hom_alt && kid.AB >= 0.98 && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq} && mom.hom_ref && mom.AB <= 0.02 && mom.DP >= ${min_dp} && mom.GQ >= ${min_gq} && dad.hom_ref && dad.AB <= 0.02 && dad.DP >= ${min_dp} && dad.GQ >= ${min_gq})" \\
    --trio "denovo_y_male:(variant.CHROM == 'chrY' && kid.sex == 1 && kid.hom_alt && kid.DP >= ${min_dp} && kid.GQ >= ${min_gq} && dad.hom_ref && dad.AB <= 0.02 && dad.DP >= ${min_dp} && dad.GQ >= ${min_gq})"

  # Add LCR annotation if the file exists
  if [ -f "${annotation_path}/LCR.bed.gz" ]; then
    echo "Adding LCR annotation"
    echo '##INFO=<ID=LCR,Number=0,Type=Flag,Description="Variant falls in low complexity region">' > lcr_header.txt
    
    # Use --mark-sites to add a flag when variants overlap with LCR regions
    bcftools annotate \
      -a ${annotation_path}/LCR.bed.gz \
      -h lcr_header.txt \
      -c CHROM,FROM,TO \
      --mark-sites +LCR \
      -O b \
      -o ${output_bcf} \
      temp_output.bcf
  else
    echo "LCR annotation file not found, skipping LCR annotation"
    # Convert output to BCF format without LCR annotation
    bcftools view -O b -o ${output_bcf} temp_output.bcf
  fi

  # Index the output BCF
  bcftools index ${output_bcf}

  # Extract all INFO field IDs from the BCF header (excluding denovo_*, highest_impact_order, CSQ, and genic)
  INFO_FIELDS=\$(bcftools view -h ${output_bcf} | grep "^##INFO=<ID=" | sed 's/##INFO=<ID=//' | cut -d',' -f1 | grep -v "^denovo_" | grep -v "^highest_impact_order\$" | grep -v "^CSQ\$" | grep -v "^genic\$" | tr '\n' ' ')

  echo "INFO fields to include in TSV: \$INFO_FIELDS"

  # Build the --info-field arguments dynamically
  INFO_ARGS=""
  for field in \$INFO_FIELDS; do
    INFO_ARGS="\$INFO_ARGS --info-field \$field"
  done

  # Extract CSQ column names from the BCF header (Format: field in Description)
  CSQ_COLUMNS=\$(bcftools view -h ${output_bcf} | grep "^##INFO=<ID=CSQ" | sed -n 's/.*Format: \\([^"]*\\).*/\\1/p' | tr '|' ' ')
  
  if [ -n "\$CSQ_COLUMNS" ]; then
    echo "CSQ columns found: \$CSQ_COLUMNS"
    # Build --csq-column arguments (one per column)
    CSQ_ARGS=""
    for column in \$CSQ_COLUMNS; do
      CSQ_ARGS="\$CSQ_ARGS --csq-column \$column"
    done
    CSQ_ARGS="\$CSQ_ARGS --csq-field CSQ"
  else
    echo "No CSQ field found in header, skipping CSQ columns"
    CSQ_ARGS=""
  fi

  # Generate TSV output from the BCF file
  slivar tsv \\
    --ped \$PED_FILE \\
    --sample-field denovo_auto_het \\
    --sample-field denovo_auto_hom \\
    --sample-field denovo_x_male \\
    --sample-field denovo_x_female_het \\
    --sample-field denovo_x_female_hom \\
    --sample-field denovo_x_par_het \\
    --sample-field denovo_x_par_hom \\
    --sample-field denovo_y_male \\
    \$INFO_ARGS \\
    \$CSQ_ARGS \\
    ${output_bcf} > ${output_tsv}
  """

  stub:
  base_name = bcf.baseName
  output_bcf = "${base_name}.dnm.bcf"
  output_tsv = "${base_name}.dnm.tsv"
  """
  touch ${output_bcf}
  touch ${output_bcf}.csi
  touch ${output_tsv}
  """
}
