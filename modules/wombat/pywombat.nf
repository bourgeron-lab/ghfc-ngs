process PYWOMBAT {
    tag "${fid}-${wombat_config_name}"
    
    publishDir "${params.data}/families/${fid}/wombat", mode: 'copy'
    
    input:
    tuple val(fid), path(annotated_parquet), path(pedigree), path(wombat_config), val(wombat_config_name), val(vep_config_name)

    output:
    tuple val(fid), val(wombat_config_name), path("${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.tsv")

    script:
    """
    wombat filter ${annotated_parquet} -p ${pedigree} -F ${wombat_config}
    """
}
