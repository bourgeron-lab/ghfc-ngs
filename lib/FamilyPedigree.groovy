/**
 * A family's own pedigree, {FID}.pedigree.tsv, held against the cohort pedigree's rows of that
 * family, and the pedigree the wombat step gives pywombat.
 *
 * FAMILIAL_PEDIGREE writes a family's file when it has none. Nothing used to write it again:
 * an edit of the cohort pedigree that keeps the family's members - a corrected father, mother or
 * sex - left the old file in place for good, and a family listed by two cohorts kept the rows of
 * whichever wrote it first. The plan now compares the two and has the file written again.
 *
 * The wombat step no longer reads the file at all. It is given the cohort's rows of the family,
 * and records them in an inputs snapshot (InputsSnapshot) beside each of its tables, so that an
 * edit of the parents or the sex rebuilds them, as it does the sampleqc and upd tables.
 *
 * Free of any Nextflow dependency, as InputsSnapshot is, so that it is tested with plain Groovy.
 */
class FamilyPedigree {

    /** pywombat's names for the columns of the pedigree it is given */
    static final String WOMBAT_HEADER = 'FID\tsample_id\tfather_id\tmother_id\tsex'

    /**
     * Each family's rows in a pedigree's lines: [fid, sample, father, mother, sex, phenotype],
     * normalised as InputsSnapshot.pedigreeRow normalises the first five, the phenotype trimmed,
     * and sorted by sample, so that a file whose rows were only reordered is the same pedigree.
     * '#' lines, blank lines and a first line starting FID are skipped. A row of fewer than six
     * columns, or a member with two different rows, is refused (IllegalArgumentException); a row
     * that repeats an earlier one of the same member is dropped.
     */
    static Map<String, List<List<String>>> rowsByFamily(List<String> lines) {
        def byFamily = [:]
        lines.eachWithIndex { String line, int index ->
            if (line.startsWith('#') || line.trim().isEmpty()) return
            def cols = line.split('\t', -1)*.trim()
            if (index == 0 && cols[0] == 'FID') return
            if (cols.size() < 6) {
                throw new IllegalArgumentException("pedigree row with ${cols.size()} columns instead of 6: ${line}")
            }
            def row = InputsSnapshot.pedigreeRow(cols) + [cols[5]]
            def members = byFamily.computeIfAbsent(row[0]) { [:] }
            def earlier = members[row[1]]
            if (earlier == null) {
                members[row[1]] = row
            } else if (earlier != row) {
                throw new IllegalArgumentException("pedigree with two different rows for ${row[1]}: ${earlier.join(' ')} and ${row.join(' ')}")
            }
        }
        return byFamily.collectEntries { fid, members ->
            [(fid): (members.values() as List).sort { a, b -> a[1] <=> b[1] }]
        }
    }

    /**
     * The rows of family fid in a family's own pedigree file, as rowsByFamily reads them, or null
     * when they cannot be read: the file is missing, empty, holds another family's rows too, or
     * is not a pedigree. A file that cannot be read is out of date.
     */
    static List<List<String>> onDisk(String path, String fid) {
        def file = new File(path)
        if (!file.isFile() || file.length() == 0) return null
        try {
            def byFamily = rowsByFamily(file.readLines('UTF-8'))
            return byFamily.keySet() == [fid] as Set ? byFamily[fid] : null
        } catch (IllegalArgumentException ignored) {
            return null
        }
    }

    /** The rows without their phenotype: what pywombat reads, and what the snapshot records. */
    static List<List<String>> withoutPhenotype(List<List<String>> rows) {
        return rows.collect { it.take(5) }
    }

    /**
     * What differs between two versions of a family's rows, one clause per member, for the run
     * summary: "KID: father DAD -> 0", "NEW: added", "OLD: removed". Empty when they are the same.
     */
    static List<String> describe(List<List<String>> before, List<List<String>> after) {
        def names = ['family', 'sample', 'father', 'mother', 'sex', 'phenotype']
        def was = (before ?: []).collectEntries { [(it[1]): it] }
        def now = (after ?: []).collectEntries { [(it[1]): it] }
        def clauses = []
        (was.keySet() + now.keySet()).sort().each { sample ->
            if (!was.containsKey(sample)) {
                clauses << "${sample}: added".toString()
            } else if (!now.containsKey(sample)) {
                clauses << "${sample}: removed".toString()
            } else {
                def changes = (2..5).findAll { was[sample][it] != now[sample][it] }
                    .collect { "${names[it]} ${was[sample][it]} -> ${now[sample][it]}" }
                if (changes) clauses << "${sample}: ${changes.join(', ')}".toString()
            }
        }
        return clauses
    }

    /**
     * The pedigree pywombat is given: a header line pywombat reads by name, then one row per
     * member without the phenotype. Unknown parents are 0, the sex 1, 2 or 0, so that a cohort
     * pedigree writing its sex as 1.0 or 2.0 reaches pywombat as 1 or 2.
     */
    static String wombatPedigree(List<List<String>> rows) {
        def text = new StringBuilder(WOMBAT_HEADER).append('\n')
        withoutPhenotype(rows).each { row -> text << row.join('\t') << '\n' }
        return text.toString()
    }
}
