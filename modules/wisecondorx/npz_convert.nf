process NPZ_CONVERT {
    tag "${barcode}"
    
    publishDir {
        def (s1, s2) = Sharding.getShards(barcode)
        "${params.data}/samples/${s1}/${s2}/${barcode}/svs/wisecondorx"
    }, mode: 'copy'
    
    input:
    tuple val(barcode), path(cram), path(crai), path(ref)
    
    output:
    tuple val(barcode), path("${barcode}.${params.wisecondorx_binsize}.npz")

    script:
    """
    wisecondorx convert ${cram} ${barcode}.${params.wisecondorx_binsize}.npz --reference ${ref} --binsize ${params.wisecondorx_binsize}
    """
}
