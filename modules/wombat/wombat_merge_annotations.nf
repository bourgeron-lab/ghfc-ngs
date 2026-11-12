process WOMBAT_MERGE_ANNOTATIONS {
    tag "${fid}-${wombat_config_name}"

    publishDir {
        def (s1, s2) = Sharding.getShards(fid)
        "${params.data}/families/${s1}/${s2}/${fid}/wombat"
    }, mode: 'copy'

    input:
    tuple val(fid), val(wombat_config_name), val(vep_config_name), path(filtered_tsv), path(vep_annotations), val(snapshot)

    output:
    tuple val(fid), val(wombat_config_name), path("${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.tsv"), emit: table
    path("inputs/${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.inputs.tsv"), emit: inputs

    script:
    // The inputs snapshot of the table, as PYWOMBAT writes it
    """
    mkdir -p inputs
    ${InputsSnapshot.decodeCommand(snapshot, "inputs/${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.inputs.tsv")}
    wombat merge-annotations \
        --tsv ${filtered_tsv} \
        --annotations ${vep_annotations} \
        -o ${fid}.rare.${vep_config_name}.annotated.${wombat_config_name} -v
    """

    stub:
    """
    mkdir -p inputs
    ${InputsSnapshot.decodeCommand(snapshot, "inputs/${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.inputs.tsv")}
    touch ${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.tsv
    """
}
