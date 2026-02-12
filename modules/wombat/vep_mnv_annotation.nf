process VEP_MNV_ANNOTATION {
    tag "${fid}-${wombat_config_name}"

    container 'ensemblorg/ensembl-vep:release_110.1'

    label 'vep'

    input:
    tuple val(fid), val(wombat_config_name), val(vep_config_name), path(vcf)

    output:
    tuple val(fid), val(wombat_config_name), val(vep_config_name), path("${fid}.mnv_annotations.${wombat_config_name}.tsv")

    script:
    """
    vep --input_file ${vcf} \
        --output_file ${fid}.mnv_annotations.${wombat_config_name}.tsv \
        --config ${params.vep_config} \
        --fork ${(1.3 * task.cpus).round()} \
        --format vcf \
        --tab \
        --force_overwrite
    """
}
