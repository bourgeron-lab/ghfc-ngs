process BCF2PARQUET {
    tag "${fid}"

    publishDir "${params.data}/families/${fid}/wombat", mode: 'copy'

    input:
    tuple val(fid), path(annotated_bcf), path(annotated_csi), val(vep_config_name)

    output:
    tuple val(fid), path("${fid}.rare.${vep_config_name}.annotated.parquet")

    script:
    """
    bcftools +split-vep -c - -p VEP_ -O b -o ${fid}.rare.${vep_config_name}.vep_splitted.bcf -d ${annotated_bcf}
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO[\\t%GT:%DP:%GQ:%AD]\\n' -HH ${fid}.rare.${vep_config_name}.vep_splitted.bcf | gzip -c > ${fid}.rare.${vep_config_name}.annotated.tsv.gz
    wombat prepare ${fid}.rare.${vep_config_name}.annotated.tsv.gz -o ${fid}.rare.${vep_config_name}.annotated.parquet
    """
}
