process PYWOMBAT {
    tag "${fid}-${wombat_config_name}"
    
    publishDir "${params.data}/families/${fid}/wombat", mode: 'copy'
    
    input:
    tuple val(fid), path(annotated_parquet), path(pedigree), path(wombat_config), val(wombat_config_name), val(vep_config_name), path(norm_bcf), path(norm_csi)

    output:
    tuple val(fid), val(wombat_config_name), path("${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.tsv")

    script:
    """
    echo "===="
    echo wombat filter ${annotated_parquet} -p ${pedigree} -F ${wombat_config} --bcf ${norm_bcf} --fasta ${params.ref}
    echo "===="
    which wombat
    echo "===="
    ls -alh /opt/conda/bin/wombat
    echo "==="
    cat ${wombat_config}
    echo "===="
    wombat filter ${annotated_parquet} -p ${pedigree} -F ${wombat_config} --bcf ${norm_bcf} --fasta ${params.ref} -v
    """
}
