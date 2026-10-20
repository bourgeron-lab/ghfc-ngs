// InputsSnapshot decides when a pedigree edit rebuilds a family's outputs: one difference in
// the text, and the family is redone; a spurious one, and every family is redone every run.

def lines = [
    'FID\tbarcode\tfather\tmother\tsex\tphenotype',
    '# a comment',
    'FAM1\tDAD\t0\t0\t1\t1',
    '',
    'FAM1\tMUM\t-9\t.\t2.0\t1',
    'FAM1\tKID\tDAD\tMUM\t1.00\t2',
    'FAM2\tSOLO\t\t0\t-9\t2 ',
]

Check.that('pedigree rows are normalised, in file order, without the phenotype') {
    assert InputsSnapshot.pedigreeRows(lines) == [
        ['FAM1', 'DAD', '0', '0', '1'],
        ['FAM1', 'MUM', '0', '0', '2'],
        ['FAM1', 'KID', 'DAD', 'MUM', '1'],
        ['FAM2', 'SOLO', '0', '0', '0'],
    ]
}
Check.that('a header is only a header on the first line') {
    assert InputsSnapshot.pedigreeRows(['FAM1\tA\t0\t0\t1\t1', 'FID\tB\t0\t0\t2\t1'])*.get(1) == ['A', 'B']
}
Check.that('a short row is refused') {
    try {
        InputsSnapshot.pedigreeRows(['FAM1\tA\t0\t0\t1'])
        assert false, 'accepted'
    }
    catch (IllegalArgumentException e) {
        assert e.message.contains('5 columns instead of 6')
    }
}
Check.that('a repeated row is dropped, even written differently or with another phenotype') {
    def repeats = lines + ['FAM1\tKID\tDAD\tMUM\t1\t1', 'FAM1\tMUM\t0\t0\t2\t-9', 'FAM1\tKID\tDAD\tMUM\t1.00\t2']
    assert InputsSnapshot.pedigreeRows(repeats) == InputsSnapshot.pedigreeRows(lines)
}
Check.that('a member with two different rows is refused') {
    ['FAM2\tKID\tDAD\tMUM\t1\t2', 'FAM1\tKID\t0\tMUM\t1\t2', 'FAM1\tKID\tDAD\tMUM\t2\t2'].each { other ->
        try {
            InputsSnapshot.pedigreeRows(lines + [other])
            assert false, "accepted ${other}"
        }
        catch (IllegalArgumentException e) {
            assert e.message.contains('two different rows for KID')
        }
    }
}
Check.that('the text is the header then the rows, one line each') {
    def text = InputsSnapshot.render([label: 'L1', mendelian: 'none'], [['F', 'A', '0', '0', '1', 'yes', '-']])
    assert text == '#inputs\tv1\n#label\tL1\n#mendelian\tnone\nF\tA\t0\t0\t1\tyes\t-\n'
}
Check.that('an edit of a parent, a sex or a member changes the text; a phenotype does not') {
    def base = InputsSnapshot.render([label: 'L'], InputsSnapshot.pedigreeRows(lines))
    def edited = { List<String> ls -> InputsSnapshot.render([label: 'L'], InputsSnapshot.pedigreeRows(ls)) }
    assert edited(lines.collect { it.replace('KID\tDAD\tMUM', 'KID\t0\tMUM') }) != base
    assert edited(lines.collect { it.replace('MUM\t-9\t.\t2.0', 'MUM\t-9\t.\t1') }) != base
    assert edited(lines + ['FAM1\tSIB\tDAD\tMUM\t2\t1']) != base
    assert edited(lines.collect { it.replace('KID\tDAD\tMUM\t1.00\t2', 'KID\tDAD\tMUM\t1\t1') }) == base
    assert InputsSnapshot.render([label: 'L2'], InputsSnapshot.pedigreeRows(lines)) != base
}
Check.that('a tab or a newline in a value is refused') {
    ['a\tb', 'a\nb'].each { bad ->
        try {
            InputsSnapshot.render([label: bad], [])
            assert false, "accepted ${bad}"
        }
        catch (IllegalArgumentException e) {
            assert e.message.contains('tab or a newline')
        }
    }
}
Check.that('the stored file matches only the same text, and a missing file never') {
    def dir = File.createTempDir()
    try {
        def text = InputsSnapshot.render([label: 'L'], [['F', 'A', '0', '0', '1']])
        def path = new File(dir, 'x.inputs.tsv').path
        assert !InputsSnapshot.matches(path, text)
        new File(path).text = text
        assert InputsSnapshot.matches(path, text)
        assert !InputsSnapshot.matches(path, text + '\n')
    }
    finally {
        dir.deleteDir()
    }
}
Check.that('the decode command writes the text byte for byte, whatever it holds') {
    def dir = File.createTempDir()
    try {
        def text = InputsSnapshot.render([label: 'L'], [['F\'1', 'A$b`c', '0', '0', '1'], ['F', 'é"\\', '0', '0', '2']])
        def path = new File(dir, 'out.tsv').path
        def proc = ['bash', '-c', InputsSnapshot.decodeCommand(text, path)].execute()
        assert proc.waitFor() == 0
        assert InputsSnapshot.matches(path, text)
    }
    finally {
        dir.deleteDir()
    }
}

Check.done()
