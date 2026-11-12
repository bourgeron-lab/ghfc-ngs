process PYWOMBAT {
    /*
    The pedigree is the cohort's rows of the family, given as text by the plan, not the family's
    {FID}.pedigree.tsv, which can be out of date. The inputs snapshot that records those rows is
    published as wombat/inputs/<table>.inputs.tsv: the plan rebuilds the table when the one it
    builds for a run differs, so an edit of a parent or a sex reruns it (lib/FamilyPedigree.groovy).
    */
    tag "${fid}-${wombat_config_name}"

    publishDir {
        def (s1, s2) = Sharding.getShards(fid)
        "${params.data}/families/${s1}/${s2}/${fid}/wombat"
    }, mode: 'copy'

    input:
    tuple val(fid), path(annotated_parquet), val(pedigree), path(wombat_config), val(wombat_config_name), val(vep_config_name), path(norm_bcf), path(norm_csi), val(snapshot)

    output:
    tuple val(fid), val(wombat_config_name), path("${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.tsv"), emit: table
    path("inputs/${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.inputs.tsv"), emit: inputs

    script:
    """
    ${InputsSnapshot.decodeCommand(pedigree, "${fid}.pedigree.tsv")}
    mkdir -p inputs
    ${InputsSnapshot.decodeCommand(snapshot, "inputs/${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.inputs.tsv")}
    echo "===="
    echo wombat filter ${annotated_parquet} -p ${fid}.pedigree.tsv -F ${wombat_config} --bcf ${norm_bcf} --fasta ${params.ref}
    echo "===="
    which wombat
    echo "===="
    ls -alh /opt/conda/bin/wombat
    echo "==="
    cat ${wombat_config}
    echo "===="
    cat ${fid}.pedigree.tsv
    echo "===="
    wombat filter ${annotated_parquet} -p ${fid}.pedigree.tsv -F ${wombat_config} --bcf ${norm_bcf} --fasta ${params.ref} -v
    """

    stub:
    """
    mkdir -p inputs
    ${InputsSnapshot.decodeCommand(snapshot, "inputs/${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.inputs.tsv")}
    touch ${fid}.rare.${vep_config_name}.annotated.${wombat_config_name}.tsv
    """
}
