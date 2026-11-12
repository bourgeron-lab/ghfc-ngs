import java.util.zip.DataFormatException
import java.util.zip.Inflater

/**
 * The samples a BCF holds, read from its header without bcftools, which is not on the launching
 * node's PATH.
 *
 * The plan needs it to tell whether a family's joint call still matches the pedigree. Inferring
 * the call's members from which of them have a gVCF only sees a member added with no gVCF: not
 * one removed from the family after the call, nor one added who already had a gVCF, nor a sample
 * replaced by another. All three leave a call that the family's later steps read as current.
 *
 * A BCF is BGZF: a run of gzip members, each with the member's size in its extra field and at most
 * 64 KiB once inflated. The header comes first - the magic `BCF\2\x`, the length of the header
 * text as a little-endian uint32, then the text - so only the first block or few are read, never
 * the records.
 *
 * Free of any Nextflow dependency, as StaleFamily is, so that it is tested with plain Groovy.
 */
class BcfHeader {

    /** No header this pipeline writes comes near it; a larger stated length is a corrupt file. */
    static final int MAX_TEXT = 64 * 1024 * 1024

    /**
     * The samples of the BCF at path, in file order, or null when it has no header to read:
     * missing, empty, truncated, not BGZF or not a BCF.
     */
    static List<String> samples(String path) {
        def file = new File(path)
        if (file.length() == 0) {
            return null
        }
        try {
            return file.withInputStream { raw ->
                def text = headerText(new DataInputStream(new BufferedInputStream(raw, 1 << 16)))
                return text == null ? null : samplesOf(text)
            }
        }
        catch (IOException | DataFormatException ignored) {
            return null
        }
    }

    /** The samples on the #CHROM line of VCF header text, or null when it has no such line. */
    static List<String> samplesOf(String text) {
        int start = 0
        if (!text.startsWith('#CHROM')) {
            int newline = text.indexOf('\n#CHROM')
            if (newline < 0) {
                return null
            }
            start = newline + 1
        }
        def end = text.indexOf('\n', start)
        // The text is NUL-terminated, and the #CHROM line is its last
        def line = (end < 0 ? text.substring(start) : text.substring(start, end)).replaceAll(/[\r\u0000]+$/, '')
        def cols = line.split('\t', -1) as List
        // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT, then one column per sample
        return cols.size() > 9 ? cols.subList(9, cols.size()) : []
    }

    // The header text, from as many blocks as it spans, or null when the file is not a BCF
    private static String headerText(DataInputStream input) {
        def out = new ByteArrayOutputStream()
        long textLength = -1
        while (true) {
            def block = nextBlock(input)
            if (block == null) {
                // The file ended before the header did
                return null
            }
            out.write(block)
            if (textLength < 0 && out.size() >= 9) {
                byte[] head = out.toByteArray()
                // B, C, F, then the major version, 2
                if (head[0] != 66 || head[1] != 67 || head[2] != 70 || head[3] != 2) {
                    return null
                }
                textLength = (head[5] & 0xffL) | ((head[6] & 0xffL) << 8) | ((head[7] & 0xffL) << 16) | ((head[8] & 0xffL) << 24)
                if (textLength > MAX_TEXT) {
                    return null
                }
            }
            if (textLength >= 0 && out.size() >= 9 + textLength) {
                return new String(out.toByteArray(), 9, (int) textLength, 'UTF-8')
            }
        }
    }

    // One BGZF block, inflated, or null at the end of the file
    private static byte[] nextBlock(DataInputStream input) {
        int first = input.read()
        if (first < 0) {
            return null
        }
        byte[] head = new byte[12]
        head[0] = (byte) first
        input.readFully(head, 1, 11)
        if ((head[0] & 0xff) != 31 || (head[1] & 0xff) != 139 || (head[2] & 0xff) != 8 || (head[3] & 4) == 0) {
            throw new IOException('not a BGZF block')
        }
        int xlen = u16(head, 10)
        byte[] extra = new byte[xlen]
        input.readFully(extra)
        // The BC subfield holds the size of the whole block, less one
        int blockSize = -1
        int i = 0
        while (i + 4 <= xlen) {
            int slen = u16(extra, i + 2)
            if ((extra[i] & 0xff) == 66 && (extra[i + 1] & 0xff) == 67 && slen == 2 && i + 6 <= xlen) {
                blockSize = u16(extra, i + 4) + 1
            }
            i += 4 + slen
        }
        // The block is the 12 bytes above, the extra field, the deflated data, a CRC and a size
        int dataLength = blockSize - 12 - xlen - 8
        if (blockSize < 0 || dataLength < 0) {
            throw new IOException('a gzip member with no BGZF block size')
        }
        byte[] deflated = new byte[dataLength]
        input.readFully(deflated)
        byte[] tail = new byte[8]
        input.readFully(tail)
        int inflatedLength = u16(tail, 4) | (u16(tail, 6) << 16)
        if (inflatedLength < 0 || inflatedLength > 65536) {
            throw new IOException('a BGZF block larger than 64 KiB')
        }
        byte[] inflated = new byte[inflatedLength]
        def inflater = new Inflater(true)
        try {
            inflater.setInput(deflated)
            int n = 0
            while (n < inflatedLength) {
                int k = inflater.inflate(inflated, n, inflatedLength - n)
                if (k == 0 && (inflater.finished() || inflater.needsInput())) {
                    break
                }
                n += k
            }
            if (n != inflatedLength) {
                throw new IOException('a BGZF block shorter than its stated size')
            }
        }
        finally {
            inflater.end()
        }
        return inflated
    }

    private static int u16(byte[] b, int at) {
        return (b[at] & 0xff) | ((b[at + 1] & 0xff) << 8)
    }
}
