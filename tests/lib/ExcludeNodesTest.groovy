// params.exclude_nodes in every form it can arrive in. nextflow.config's process.ext.slurm_options
// splits the value the same way for the tasks' --exclude=; tests/stub/run.sh checks that side.

Check.that('none: empty, missing, false') {
    assert ExcludeNodes.parse('') == []
    assert ExcludeNodes.parse(null) == []
    assert ExcludeNodes.parse(false) == []
    assert ExcludeNodes.parse([]) == []
    ['', null, false, []].each { assert ExcludeNodes.check(it) == null }
}
Check.that('a comma-separated string, as from the command line') {
    assert ExcludeNodes.parse('maestro-1096,maestro-1082') == ['maestro-1096', 'maestro-1082']
    assert ExcludeNodes.parse('maestro-1096') == ['maestro-1096']
    assert ExcludeNodes.check('maestro-1096,maestro-1082') == null
}
Check.that('a YAML list, a number, stray commas and spaces') {
    assert ExcludeNodes.parse(['maestro-1096', 'maestro-1082']) == ['maestro-1096', 'maestro-1082']
    assert ExcludeNodes.parse(1096) == ['1096']
    assert ExcludeNodes.parse(' maestro-1096, ,maestro-1082,') == ['maestro-1096', 'maestro-1082']
    assert ExcludeNodes.check(['maestro-1096', 'node.a_b']) == null
}
Check.that('a bare flag is refused, not read as no node') {
    assert ExcludeNodes.parse(true) == []
    assert ExcludeNodes.check(true).contains('no value')
}
Check.that('names Slurm would not take are refused, by name') {
    assert ExcludeNodes.check('maestro-1096,bad node').contains("'bad node' is not a node name")
    assert ExcludeNodes.check('maestro-[1082,1096]').contains('are not node names')
    assert ExcludeNodes.check('a;rm -rf x') != null
    assert ExcludeNodes.check([['a']]) != null
    assert ExcludeNodes.check([a: 1]) != null
}

Check.done()
