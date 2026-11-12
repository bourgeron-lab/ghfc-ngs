// BcfHeader reads which samples a family's norm.bcf holds, for the plan's drift check. The
// fixtures were written by bcftools 1.22 (view --no-version), so the reader is checked against
// real BGZF and not only against itself:
//   trio.bcf         KID, FATHER, MOTHER, one record
//   long_header.bcf  SOLO, a ~300 KB header of 4000 contigs, over several BGZF blocks
//   sites_only.bcf   no samples
//   trio.vcf.gz      the trio as a bgzipped VCF, which is BGZF but not a BCF

def fixtures = new File(new File(getClass().protectionDomain.codeSource.location.toURI()).parentFile, 'fixtures')
def fixture = { String name -> new File(fixtures, name).path }

Check.that('a BCF gives its samples in file order') {
    assert BcfHeader.samples(fixture('trio.bcf')) == ['KID', 'FATHER', 'MOTHER']
}
Check.that('a header over several BGZF blocks is read whole') {
    assert BcfHeader.samples(fixture('long_header.bcf')) == ['SOLO']
}
Check.that('a BCF with no sample columns has no samples, which is not unreadable') {
    assert BcfHeader.samples(fixture('sites_only.bcf')) == []
}
Check.that('a file that is not a BCF has no header to read') {
    assert BcfHeader.samples(fixture('trio.vcf.gz')) == null
    def dir = File.createTempDir()
    try {
        def text = new File(dir, 'text.bcf')
        text.text = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tKID\n'
        assert BcfHeader.samples(text.path) == null
    }
    finally { dir.deleteDir() }
}
Check.that('a missing, empty or truncated BCF has no header to read') {
    def dir = File.createTempDir()
    try {
        assert BcfHeader.samples(new File(dir, 'absent.bcf').path) == null
        def empty = new File(dir, 'empty.bcf')
        empty.bytes = new byte[0]
        assert BcfHeader.samples(empty.path) == null
        // Cut inside the header's blocks, and inside the very first one
        def whole = new File(fixture('long_header.bcf')).bytes
        [whole.length.intdiv(4), 40, 5].each { n ->
            def cut = new File(dir, "cut${n}.bcf")
            cut.bytes = Arrays.copyOf(whole, n)
            assert BcfHeader.samples(cut.path) == null
        }
    }
    finally { dir.deleteDir() }
}
Check.that('the samples come from the #CHROM line, wherever it is and however it ends') {
    def chrom = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
    assert BcfHeader.samplesOf("##fileformat=VCFv4.2\n${chrom}\tA\tB\n\u0000") == ['A', 'B']
    assert BcfHeader.samplesOf("##fileformat=VCFv4.2\n${chrom}\tA\tB\u0000") == ['A', 'B']
    assert BcfHeader.samplesOf("${chrom}\tA\r\n") == ['A']
    assert BcfHeader.samplesOf("##x\n${chrom}\n") == []
    assert BcfHeader.samplesOf("##x\n#CHROM\tPOS\n") == []
    assert BcfHeader.samplesOf("##fileformat=VCFv4.2\n##comment #CHROM\n") == null
}

Check.done()
