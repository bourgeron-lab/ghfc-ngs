process PYWOMBAT {
    tag "${fid}-${wombat_config_name}"
    
    publishDir "${params.data}/families/${fid}/wombat", mode: 'copy'
    
    input:
    tuple val(fid), path(annotated_tsv), path(pedigree), path(wombat_config), val(wombat_config_name), val(vep_config_name)
    
    output:
    tuple val(fid), val(wombat_config_name), path("${fid}.rare.${vep_config_name}.${wombat_config_name}.tsv.gz")
    
    script:
    """
    wombat ${annotated_tsv} -p ${pedigree} -F ${wombat_config}
    """
}
