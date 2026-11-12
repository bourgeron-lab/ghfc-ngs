// FamilyPedigree decides when a family's own pedigree file is out of date, and what pywombat is
// given: a difference it misses leaves a family analysed with the wrong parents for good; one it
// invents rewrites every family's file, and reruns every family's wombat tables, every run.

def cohort = [
    'FAM1\tKID\tDAD\tMUM\t1\t2',
    'FAM1\tDAD\t0\t0\t1\t1',
    'FAM1\tMUM\t0\t0\t2\t1',
    'FAM2\tSOLO\t0\t0\t2.0\t2',
]

Check.that('rows are grouped by family, normalised and sorted by sample') {
    def rows = FamilyPedigree.rowsByFamily(cohort)
    assert rows.keySet() == ['FAM1', 'FAM2'] as Set
    assert rows.FAM1 == [
        ['FAM1', 'DAD', '0', '0', '1', '1'],
        ['FAM1', 'KID', 'DAD', 'MUM', '1', '2'],
        ['FAM1', 'MUM', '0', '0', '2', '1'],
    ]
    assert rows.FAM2 == [['FAM2', 'SOLO', '0', '0', '2', '2']]
}
Check.that('a file whose rows are only reordered, or written differently, is the same pedigree') {
    def other = ['FID\tIID\tPAT\tMAT\tSEX\tPHENO', '# note', 'FAM1\tMUM\t-9\t.\t2.0\t1', '',
                 'FAM1\tKID\tDAD\tMUM\t1.00\t2 ', 'FAM1\tDAD\t\t0\t1\t1', 'FAM1\tDAD\t0\t0\t1\t1']
    assert FamilyPedigree.rowsByFamily(other).FAM1 == FamilyPedigree.rowsByFamily(cohort).FAM1
}
Check.that('an edit of a parent, a sex or a phenotype is a difference') {
    def base = FamilyPedigree.rowsByFamily(cohort).FAM1
    ['FAM1\tKID\tDAD\t0\t1\t2', 'FAM1\tKID\tDAD\tMUM\t2\t2', 'FAM1\tKID\tDAD\tMUM\t1\t1'].each { edit ->
        def edited = FamilyPedigree.rowsByFamily([edit] + cohort.drop(1)).FAM1
        assert edited != base
    }
}
Check.that('a short row, or a member with two different rows, is refused') {
    [['FAM1\tA\t0\t0\t1'], ['FAM1\tA\t0\t0\t1\t1', 'FAM1\tA\t0\t0\t2\t1']].each { lines ->
        try {
            FamilyPedigree.rowsByFamily(lines)
            assert false, "accepted ${lines}"
        } catch (IllegalArgumentException ignored) {
        }
    }
}

def dir = File.createTempDir()
def write = { String name, List<String> lines ->
    def f = new File(dir, name)
    f.text = lines.collect { it + '\n' }.join('')
    f.path
}

Check.that('a family file is read as its own family\'s rows, header or not') {
    assert FamilyPedigree.onDisk(write('a.tsv', cohort.take(3)), 'FAM1') == FamilyPedigree.rowsByFamily(cohort).FAM1
    assert FamilyPedigree.onDisk(write('b.tsv', ['FID\tIID\tPAT\tMAT\tSEX\tPHENO'] + cohort.take(3)), 'FAM1') ==
        FamilyPedigree.rowsByFamily(cohort).FAM1
}
Check.that('a family file that cannot be read as the family\'s pedigree is null') {
    assert FamilyPedigree.onDisk(new File(dir, 'absent.tsv').path, 'FAM1') == null
    assert FamilyPedigree.onDisk(write('empty.tsv', []), 'FAM1') == null
    assert FamilyPedigree.onDisk(write('short.tsv', ['FAM1\tKID\tDAD\tMUM\t1']), 'FAM1') == null
    assert FamilyPedigree.onDisk(write('other.tsv', cohort), 'FAM1') == null
    assert FamilyPedigree.onDisk(write('wrong.tsv', cohort.drop(3)), 'FAM1') == null
}
Check.that('describe says what each member gained, lost or changed') {
    def before = FamilyPedigree.rowsByFamily([
        'F\tA\t0\t0\t1\t1', 'F\tB\tX\tY\t2\t1', 'F\tGONE\t0\t0\t1\t1', 'F\tSAME\t0\t0\t2\t1']).F
    def after = FamilyPedigree.rowsByFamily([
        'F\tA\t0\t0\t1\t2', 'F\tB\tX\t0\t1\t1', 'F\tNEW\t0\t0\t1\t1', 'F\tSAME\t0\t0\t2\t1']).F
    assert FamilyPedigree.describe(before, after) == [
        'A: phenotype 1 -> 2',
        'B: mother Y -> 0, sex 2 -> 1',
        'GONE: removed',
        'NEW: added',
    ]
    assert FamilyPedigree.describe(after, after) == []
}
Check.that('pywombat gets a header it reads by name, then the rows without phenotype, sex as 1 or 2') {
    assert FamilyPedigree.wombatPedigree(FamilyPedigree.rowsByFamily(cohort).FAM2) ==
        'FID\tsample_id\tfather_id\tmother_id\tsex\nFAM2\tSOLO\t0\t0\t2\n'
    assert FamilyPedigree.wombatPedigree(FamilyPedigree.rowsByFamily(cohort).FAM1).readLines()[1] == 'FAM1\tDAD\t0\t0\t1'
}
Check.that('the wombat snapshot ignores the phenotype, so a phenotype edit rebuilds no table') {
    def rows = FamilyPedigree.rowsByFamily(cohort).FAM1
    def edited = FamilyPedigree.rowsByFamily(['FAM1\tKID\tDAD\tMUM\t1\t1'] + cohort.drop(1)).FAM1
    assert FamilyPedigree.withoutPhenotype(rows) == FamilyPedigree.withoutPhenotype(edited)
}

dir.deleteDir()

Check.done()
