/**
 * A family's own pedigree, {FID}.pedigree.tsv, held against the cohort pedigree's rows of that
 * family, and the pedigree the wombat step gives pywombat.
 *
 * FAMILIAL_PEDIGREE writes a family's file when it has none. Nothing used to write it again:
 * an edit of the cohort pedigree that keeps the family's members - a corrected father, mother or
 * sex - left the old file in place for good, and a family listed by two cohorts kept the rows of
 * whichever wrote it first. The plan now compares the two and has the file written again.
 *
 * Unless another cohort of the data directory lists the family with other rows (conflicts):
 * then each cohort's run would write the file, and redo the family's outputs, in turn. Such a
 * family is held by default (main.nf, pedigree_conflicts), and its file left as it is.
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

    // ------------------------------------------------------------------------------------
    // The other cohorts of a data directory
    // ------------------------------------------------------------------------------------
    //
    // The cohorts of a data directory share one families/ tree, but each has its own pedigree.
    // A family two cohorts list with different rows would be written by each cohort's run in
    // turn, {FID}.pedigree.tsv and every family-level output, and belong to whichever ran last.

    /**
     * The pedigrees of the data directory's other cohorts: every cohorts/<NAME>/<NAME>.pedigree.tsv
     * but this cohort's, and but the file this cohort's run reads (ownPedigree, may be null).
     * Hidden directories are skipped. Returns [rows: cohort -> rowsByFamily, unreadable: cohort ->
     * why], each sorted by cohort: a pedigree that cannot be read is reported, not fatal.
     */
    static Map otherCohorts(String data, String cohort, String ownPedigree = null) {
        def rows = new TreeMap<String, Map<String, List<List<String>>>>()
        def unreadable = new TreeMap<String, String>()
        def dirs = new File(data, 'cohorts').listFiles()?.findAll { it.isDirectory() && !it.name.startsWith('.') } ?: []
        def own = ownPedigree ? canonical(new File(ownPedigree)) : null
        dirs.each { File dir ->
            if (dir.name == cohort) return
            def file = new File(dir, "${dir.name}.pedigree.tsv")
            if (!file.isFile() || (own != null && canonical(file) == own)) return
            try {
                rows[dir.name] = rowsByFamily(file.readLines('UTF-8'))
            } catch (IllegalArgumentException | IOException e) {
                unreadable[dir.name] = "${file.path}: ${e.message}".toString()
            }
        }
        return [rows: rows, unreadable: unreadable]
    }

    private static String canonical(File file) {
        try {
            return file.canonicalPath
        } catch (IOException ignored) {
            return file.absolutePath
        }
    }

    /**
     * The families of this cohort (mine, as rowsByFamily) that another cohort lists with other
     * rows once the phenotype is dropped: another father, mother, sex, member or dummy. A family
     * whose rows differ only in the phenotype is not one: two cohorts may study one family for two
     * different conditions. One entry per family, sorted:
     *
     *   [family: 'FAM01', cohorts: ['OTHER'], changes: [OTHER: ['KID: father 0 -> DAD']]]
     *
     * Each clause, as describe words it, reads this cohort's rows -> the other cohort's.
     */
    static List<Map> conflicts(Map<String, List<List<String>>> mine, Map<String, Map<String, List<List<String>>>> others) {
        def found = []
        (mine.keySet() as List).sort().each { fid ->
            def ours = withoutPhenotype(mine[fid])
            def changes = new TreeMap<String, List<String>>()
            others.each { cohort, families ->
                def theirs = families[fid]
                if (theirs == null) return
                theirs = withoutPhenotype(theirs)
                if (theirs != ours) changes[cohort] = describe(ours, theirs)
            }
            if (changes) {
                found << [family: fid, cohorts: changes.keySet() as List, changes: changes as LinkedHashMap]
            }
        }
        return found
    }

    /** How the run summary words a conflict: "FAM01 differs from OTHER: KID: father 0 -> DAD" */
    static List<String> describeConflict(Map conflict) {
        return conflict.changes.collect { cohort, clauses ->
            "${conflict.family} differs from ${cohort}: ${clauses.join('; ')}".toString()
        }
    }

    /**
     * The families of this cohort whose samples another cohort lists under another family ID:
     * the same people twice, so their family-level work is done twice. A dummy is not a sample,
     * so one shared by two families is not an alias. One entry per family and other cohort's
     * family, sorted: [family: 'FAM02', cohort: 'OTHER', as: 'FAM99', samples: ['BC004']].
     */
    static List<Map> aliases(Map<String, List<List<String>>> mine, Map<String, Map<String, List<List<String>>>> others) {
        // Each other cohort's sample -> its family, once
        def familyOf = new TreeMap<String, Map<String, String>>()
        others.each { cohort, families ->
            def index = [:]
            families.each { otherFid, rows ->
                rows.each { row -> if (!DummyMember.isDummy(row[1])) index[row[1]] = otherFid }
            }
            familyOf[cohort] = index
        }
        def found = []
        (mine.keySet() as List).sort().each { fid ->
            def samples = mine[fid].collect { it[1] }.findAll { !DummyMember.isDummy(it) }
            familyOf.each { cohort, index ->
                def shared = new TreeMap<String, List<String>>()
                samples.each { sample ->
                    def otherFid = index[sample]
                    if (otherFid != null && otherFid != fid) shared.computeIfAbsent(otherFid) { [] } << sample
                }
                shared.each { otherFid, list ->
                    found << [family: fid, cohort: cohort, as: otherFid, samples: list.sort()]
                }
            }
        }
        return found
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
