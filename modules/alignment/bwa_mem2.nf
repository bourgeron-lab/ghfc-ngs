process BWA_MEM2_ALIGN {
    
    tag "$barcode:$unit"
    
    publishDir "${params.scratch}/cram", mode: 'copy'
    
    input:
    tuple val(barcode), val(unit), path(fastq1), path(fastq2), val(project), val(flowcell), val(dual), val(lane)
    val ref
    val bwa_mem2_path
    val samblaster_path
    val sambamba_path
    val samtools_path
    val scratch_dir
    
    output:
    tuple val(barcode), val(unit), path("${unit}.cram"), emit: cram
    
    script:
    def read_group = "@RG\\tID:A${project}.DA.${barcode}.${flowcell}.${dual}.${lane}\\tDS:A${project}\\tSM:${barcode}\\tPL:ILLUMINA"
    
    """
    mkdir -p ${scratch_dir}/tmp
    
    ${bwa_mem2_path} mem -Y -t ${task.cpus - 4} \\
        -R "${read_group}" \\
        ${ref} \\
        ${fastq1} \\
        ${fastq2} \\
        | ${samblaster_path} --addMateTags \\
        | ${sambamba_path} view -S -f bam -l 0 -t 4 -o /dev/stdout /dev/stdin \\
        | ${sambamba_path} sort --tmpdir=${scratch_dir}/tmp/ -l 0 -m 50GiB -t ${task.cpus} -o /dev/stdout /dev/stdin \\
        | ${samtools_path} view -C -@${task.cpus} -T ${ref} -o ${unit}.cram -
    """
    
    stub:
    """
    touch ${unit}.cram
    """
}
