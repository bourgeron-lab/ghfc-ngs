process PYWOMBAT_FILTER {
    tag "${fid}-${wombat_config_name}"

    input:
    tuple val(fid), path(annotated_parquet), path(pedigree), path(wombat_config), val(wombat_config_name), val(vep_config_name), path(norm_bcf), path(norm_csi)

    output:
    tuple val(fid), val(wombat_config_name), val(vep_config_name), path("${fid}.rare.${vep_config_name}.filtered.${wombat_config_name}.tsv"), path("${fid}.rare.${vep_config_name}.filtered.${wombat_config_name}.to_annotate.vcf")

    script:
    """
    wombat filter ${annotated_parquet} -p ${pedigree} -F ${wombat_config} \
        --bcf ${norm_bcf} --fasta ${params.ref} \
        -o ${fid}.rare.${vep_config_name}.filtered.${wombat_config_name} -v
    """
}
