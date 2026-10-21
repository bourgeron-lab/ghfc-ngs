/**
 * A record of the inputs a family's (or the cohort's) outputs were built from, written beside
 * them, so that the plan can tell a pedigree edit from a tree that is merely complete.
 *
 * Every other step decides what to rebuild from which files exist. That cannot see an edit
 * that keeps a family's members: a corrected father, mother or sex leaves every output in
 * place and wrong. A step that uses the pedigree that way publishes the text built here as
 * <prefix>.inputs.tsv, and the plan rebuilds whenever the text it builds for this run is not
 * byte for byte the one on disk.
 *
 * The text is built here at planning time and handed to the task, which writes it verbatim
 * (decodeCommand): the stored and the expected snapshot are made by the same code, so a
 * normalisation in a shell script can never make every family look changed.
 *
 * The format is TSV: '#' lines are a header of key and value, then one row per member. The
 * reading scripts skip the header, so it can carry anything that should force a rebuild when
 * it changes.
 *
 * Free of any Nextflow dependency, as StaleFamily is, so that it is tested with plain Groovy.
 */
class InputsSnapshot {

    static final String VERSION = 'v1'

    /**
     * A pedigree file's rows, normalised and in file order: [fid, sample, father, mother, sex].
     * Read as main.nf's parsePedigreeFile reads it: '#' lines and blank lines skipped, and a
     * first line starting FID taken for a header. An unknown parent (0, -9, . or empty) is 0;
     * the sex is 1 or 2, with a trailing .0 dropped, and 0 for anything else. The phenotype is
     * left out: no step that uses this reads it, and an edit of it should rebuild nothing.
     * A row that repeats an earlier one of the same member is dropped, as parsePedigreeFile
     * drops it; a member given two different rows is refused.
     */
    static List<List<String>> pedigreeRows(List<String> lines) {
        def rows = [:]
        lines.eachWithIndex { String line, int index ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def cols = line.split('\t', -1)*.trim()
            if (index == 0 && cols[0] == 'FID') return
            if (cols.size() < 6) {
                throw new IllegalArgumentException("pedigree row with ${cols.size()} columns instead of 6: ${line}")
            }
            def row = pedigreeRow(cols)
            def earlier = rows[row[1]]
            if (earlier == null) {
                rows[row[1]] = row
            } else if (earlier != row) {
                throw new IllegalArgumentException("pedigree with two different rows for ${row[1]}: ${earlier.join(' ')} and ${row.join(' ')}")
            }
        }
        return rows.values() as List
    }

    /**
     * One pedigree row's first five columns, normalised as pedigreeRows says: two rows of a
     * member that only differ in how they write an unknown parent or a sex are the same row.
     */
    static List<String> pedigreeRow(List<String> cols) {
        def parent = { String p -> p in ['', '0', '-9', '.'] ? '0' : p }
        def sex = cols[4].trim().replaceAll(/\.0+$/, '')
        return [cols[0].trim(), cols[1].trim(), parent(cols[2].trim()), parent(cols[3].trim()), sex in ['1', '2'] ? sex : '0']
    }

    /**
     * The snapshot text: the header entries in their order, then the rows, every line ending
     * in a newline. A tab or a newline in a value would make another file, so it is refused.
     */
    static String render(Map<String, String> header, List<List<String>> rows) {
        def check = { String value ->
            if (value == null || value.contains('\t') || value.contains('\n') || value.contains('\r')) {
                throw new IllegalArgumentException("a snapshot value may not be null, or hold a tab or a newline: ${value}")
            }
            return value
        }
        def text = new StringBuilder("#inputs\t${VERSION}\n")
        header.each { key, value -> text << "#${check(key)}\t${check(value)}\n" }
        rows.each { row -> text << row.collect { check(it?.toString()) }.join('\t') << '\n' }
        return text.toString()
    }

    /** Whether the file holds exactly this snapshot. A missing file never does. */
    static boolean matches(String path, String snapshot) {
        def stored = new File(path)
        return stored.isFile() && stored.getText('UTF-8') == snapshot
    }

    /**
     * A shell command that writes the snapshot to a file, byte for byte. Base64 keeps every
     * character of a sample name out of the shell's reach.
     */
    static String decodeCommand(String snapshot, String path) {
        def encoded = snapshot.getBytes('UTF-8').encodeBase64().toString()
        return "printf '%s' '${encoded}' | base64 -d > '${path}'"
    }
}
