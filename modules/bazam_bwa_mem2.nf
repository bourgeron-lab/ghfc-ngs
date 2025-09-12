process BAZAM_BWA_MEM2_REALIGN {

    tag "$barcode"
    
    publishDir "${data_dir}/cram", mode: 'copy'

    input:
    tuple val(barcode), path(input_cram), path(input_crai)
    val ref
    val oldref
    val bwa_mem2_path
    val bazam_path
    val samblaster_path
    val sambamba_path
    val samtools_path
    val scratch_dir
    val data_dir
    val ref_name

    output:
    tuple val(barcode), path("${barcode}.${ref_name}.cram"), path("${barcode}.${ref_name}.cram.crai"), emit: cram

    script:
    def read_group = "@RG\\tID:${barcode}\\tSM:${barcode}\\tPL:ILLUMINA"

    """
    mkdir -p ${scratch_dir}/tmp
    set -euo pipefail

    java -Xmx12g -Dsamjdk.reference_fasta=${oldref} -jar ${bazam_path} -bam ${input_cram} \\
        | ${bwa_mem2_path} mem -Y -t ${task.cpus - 4} -p \\
            -R "${read_group}" \\
            ${ref} - \\
        | ${samblaster_path} --addMateTags --ignoreUnmated \\
        | ${sambamba_path} view -S -f bam -l 0 -t 4 -o /dev/stdout /dev/stdin \\
        | ${sambamba_path} sort --tmpdir=${scratch_dir}/tmp/ -l 0 -m 50GiB -t ${task.cpus} -o /dev/stdout /dev/stdin \\
        | ${samtools_path} view -C -@${task.cpus} -T ${ref} -o ${barcode}.${ref_name}.cram -

    ${samtools_path} index ${barcode}.${ref_name}.cram
    """

    stub:
    """
    touch ${barcode}.${ref_name}.cram
    touch ${barcode}.${ref_name}.cram.crai
    """
}
