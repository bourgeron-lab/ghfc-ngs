process NORMALIZE {

    tag "$fid"

    publishDir "${params.data}/families/${fid}/vcfs", mode: 'copy'

    input:
    tuple val(fid), path(vcf), path(tbi)

    output:
    tuple val(fid), path("${fid}.norm.vcf.gz"), path("${fid}.norm.vcf.gz.tbi"), emit: normalized_vcf

    script:
    """
    # Normalize the VCF file
    zcat ${vcf} | \\
        sed -e 's/ID=AD,Number=\\./ID=AD,Number=R/' | \\
        bcftools norm -m - -w 1000 -c s -f ${params.ref} --threads ${task.cpus} -Oz -o ${fid}.norm.vcf.gz

    # Index the normalized VCF
    bcftools index --tbi ${fid}.norm.vcf.gz
    """

    stub:
    """
    touch ${fid}.norm.vcf.gz
    touch ${fid}.norm.vcf.gz.tbi
    """
}