class Sharding {

    static List getShards(String id) {
        def stripped = id.replaceAll(/[-._]/, '')
        def shard1 = stripped[-1].toUpperCase()
        def shard2 = stripped[-2].toUpperCase()
        return [shard1, shard2]
    }

    static String getSampleDir(String data, String barcode) {
        def (s1, s2) = getShards(barcode)
        return "${data}/samples/${s1}/${s2}/${barcode}"
    }

    static String getFamilyDir(String data, String fid) {
        def (s1, s2) = getShards(fid)
        return "${data}/families/${s1}/${s2}/${fid}"
    }
}
