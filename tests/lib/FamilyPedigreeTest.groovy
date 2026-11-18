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

// The other cohorts of a data directory: a family this cohort lists that another lists with
// other rows is held rather than overwritten (main.nf, pedigree_conflicts). A difference missed
// lets two cohorts rewrite one family's outputs in turn; one invented holds a family for nothing.

def mine = FamilyPedigree.rowsByFamily([
    'FAM01\tKID\t0\tMUM\t1\t2', 'FAM01\tDAD\t0\t0\t1\t1', 'FAM01\tMUM\t0\t0\t2\t1',
    'FAM02\tSOLO\t0\t0\t2\t2',
    'FAM03\tALONE\t0\t0\t1\t2',
    'FAM04\tTWIN\tDUMMY_TWIN_F\tDUMMY_TWIN_M\t1\t2', 'FAM04\tDUMMY_TWIN_F\t0\t0\t1\t0', 'FAM04\tDUMMY_TWIN_M\t0\t0\t2\t0',
])

Check.that('a family another cohort lists with other parents, sex or members is a conflict, worded this cohort -> the other') {
    def others = [
        OTHER: FamilyPedigree.rowsByFamily([
            'FAM01\tKID\tDAD\tMUM\t1\t2', 'FAM01\tDAD\t0\t0\t1\t1', 'FAM01\tMUM\t0\t0\t2\t1',
            'FAM02\tSOLO\t0\t0\t1\t2']),
        ALSO: FamilyPedigree.rowsByFamily([
            'FAM01\tKID\t0\tMUM\t1\t2', 'FAM01\tMUM\t0\t0\t2\t1']),
    ]
    def found = FamilyPedigree.conflicts(mine, others)
    assert found == [
        [family: 'FAM01', cohorts: ['ALSO', 'OTHER'],
         changes: [ALSO: ['DAD: removed'], OTHER: ['KID: father 0 -> DAD']]],
        [family: 'FAM02', cohorts: ['OTHER'], changes: [OTHER: ['SOLO: sex 2 -> 1']]],
    ]
    assert FamilyPedigree.describeConflict(found[0]) == ['FAM01 differs from ALSO: DAD: removed',
                                                         'FAM01 differs from OTHER: KID: father 0 -> DAD']
}
Check.that('a phenotype, the row order or the way a row is written is not a conflict') {
    def others = [OTHER: FamilyPedigree.rowsByFamily([
        'FID\tIID\tPAT\tMAT\tSEX\tPHENO',
        'FAM01\tMUM\t-9\t.\t2.0\t2', 'FAM01\tKID\t0\tMUM\t1\t1', 'FAM01\tDAD\t\t0\t1\t1',
        'FAM03\tALONE\t0\t0\t1\t1'])]
    assert FamilyPedigree.conflicts(mine, others) == []
}
Check.that('dummies count: a family whose dummy parents another cohort lacks or names otherwise is a conflict') {
    def others = [
        BARE: FamilyPedigree.rowsByFamily(['FAM04\tTWIN\t0\t0\t1\t2']),
        RENAMED: FamilyPedigree.rowsByFamily([
            'FAM04\tTWIN\tDUMMY_TWIN_F\tDUMMY_TWIN_MOTHER\t1\t2', 'FAM04\tDUMMY_TWIN_F\t0\t0\t1\t0',
            'FAM04\tDUMMY_TWIN_MOTHER\t0\t0\t2\t0']),
    ]
    assert FamilyPedigree.conflicts(mine, others) == [
        [family: 'FAM04', cohorts: ['BARE', 'RENAMED'],
         changes: [BARE: ['DUMMY_TWIN_F: removed', 'DUMMY_TWIN_M: removed', 'TWIN: father DUMMY_TWIN_F -> 0, mother DUMMY_TWIN_M -> 0'],
                   RENAMED: ['DUMMY_TWIN_M: removed', 'DUMMY_TWIN_MOTHER: added', 'TWIN: mother DUMMY_TWIN_M -> DUMMY_TWIN_MOTHER']]],
    ]
}
Check.that('a family no other cohort lists, or no other cohort at all, is no conflict') {
    assert FamilyPedigree.conflicts(mine, [:]) == []
    assert FamilyPedigree.conflicts(mine, [OTHER: FamilyPedigree.rowsByFamily(['FAM99\tX\t0\t0\t1\t1'])]) == []
}

Check.that('a sample another cohort lists under another family ID is an alias; a shared dummy is not') {
    def others = [
        OTHER: FamilyPedigree.rowsByFamily([
            'FAM99\tKID\tDAD\tMUM\t1\t2', 'FAM99\tDAD\t0\t0\t1\t1', 'FAM99\tMUM\t0\t0\t2\t1',
            'FAM98\tSOLO\t0\t0\t2\t2',
            'FAM97\tOTHERTWIN\tDUMMY_TWIN_F\tDUMMY_TWIN_M\t1\t2', 'FAM97\tDUMMY_TWIN_F\t0\t0\t1\t0', 'FAM97\tDUMMY_TWIN_M\t0\t0\t2\t0',
            'FAM03\tALONE\t0\t0\t1\t2']),
        ALSO: FamilyPedigree.rowsByFamily(['FAM01\tKID\t0\tMUM\t1\t2', 'X\tMUM\t0\t0\t2\t1']),
    ]
    assert FamilyPedigree.aliases(mine, others) == [
        [family: 'FAM01', cohort: 'ALSO', as: 'X', samples: ['MUM']],
        [family: 'FAM01', cohort: 'OTHER', as: 'FAM99', samples: ['DAD', 'KID', 'MUM']],
        [family: 'FAM02', cohort: 'OTHER', as: 'FAM98', samples: ['SOLO']],
    ]
}

def data = File.createTempDir()
def cohortFile = { String name, List<String> lines ->
    def d = new File(data, "cohorts/${name}")
    d.mkdirs()
    def f = new File(d, "${name}.pedigree.tsv")
    f.text = lines.collect { it + '\n' }.join('')
    f.path
}

Check.that('the other cohorts are every cohorts/<NAME>/<NAME>.pedigree.tsv but this one and hidden ones; one that cannot be read is reported') {
    def own = cohortFile('MINE', ['FAM01\tKID\t0\t0\t1\t2'])
    cohortFile('OTHER', ['FAM01\tKID\tDAD\t0\t1\t2'])
    cohortFile('.hidden', ['FAM01\tKID\tX\t0\t1\t2'])
    cohortFile('BROKEN', ['FAM01\tKID\t0\t0\t1'])
    new File(data, 'cohorts/EMPTYDIR').mkdirs()
    new File(data, 'cohorts/WRONGNAME').mkdirs()
    new File(data, 'cohorts/WRONGNAME/other.pedigree.tsv').text = 'FAM01\tKID\tY\t0\t1\t2\n'
    def found = FamilyPedigree.otherCohorts(data.path, 'MINE', own)
    assert found.rows.keySet() as List == ['OTHER']
    assert found.rows.OTHER.FAM01 == [['FAM01', 'KID', 'DAD', '0', '1', '2']]
    assert found.unreadable.keySet() as List == ['BROKEN']
    assert found.unreadable.BROKEN.contains('5 columns instead of 6')
}
Check.that('a pedigree kept in another cohort\'s directory is not compared with itself') {
    def borrowed = new File(data, 'cohorts/OTHER/OTHER.pedigree.tsv').path
    assert !('OTHER' in FamilyPedigree.otherCohorts(data.path, 'ELSEWHERE', borrowed).rows.keySet())
    assert FamilyPedigree.otherCohorts(new File(data, 'absent').path, 'MINE').rows == [:]
}

data.deleteDir()

Check.done()
