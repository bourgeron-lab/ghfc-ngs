process NPZ_CONVERT {
    tag "${barcode}"
    
    publishDir {
        def (s1, s2) = Sharding.getShards(barcode)
        "${params.data}/samples/${s1}/${s2}/${barcode}/svs/wisecondorx"
    }, mode: 'copy'
    
    input:
    tuple val(barcode), path(cram), path(crai), path(ref)
    
    output:
    tuple val(barcode), path("${barcode}.npz")
    
    script:
    """
    WisecondorX convert ${cram} ${barcode}.npz --reference ${ref}
    """
}
