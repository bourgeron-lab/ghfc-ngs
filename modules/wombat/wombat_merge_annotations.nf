process WOMBAT_MERGE_ANNOTATIONS {
    tag "${fid}-${wombat_config_name}"

    publishDir {
        def (s1, s2) = Sharding.getShards(fid)
        "${params.data}/families/${s1}/${s2}/${fid}/wombat"
    }, mode: 'copy'

    input:
    tuple val(fid), val(wombat_config_name), val(vep_config_name), path(filtered_tsv), path(vep_annotations)

    output:
    tuple val(fid), val(wombat_config_name), path("${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.tsv")

    script:
    """
    wombat merge-annotations \
        --tsv ${filtered_tsv} \
        --annotations ${vep_annotations} \
        -o ${fid}.rare.${vep_config_name}.annotated.${wombat_config_name} -v
    """
}
