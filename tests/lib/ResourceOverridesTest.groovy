// The startup check of a cohort's resources: block. A mistake it lets through surfaces hours
// into a run, or never.

def root = new File(getClass().protectionDomain.codeSource.location.path).parentFile.parentFile.parentFile
def known = ResourceOverrides.processNames(root)
def ceilings = [memory: '500.GB', cpus: 96, time: '48.h']
def check = { Map resources -> ResourceOverrides.check(resources, known, ceilings) }

Check.that('processNames finds the modules and the include aliases') {
    assert known.names.PYWOMBAT == 'PYWOMBAT'
    assert known.names.ROH_CALL == 'ROH_CALL'
    assert known.aliases.values().every { it in known.names.keySet() }
}
Check.that('nothing to check is fine') {
    assert check(null) == [errors: [], warnings: []]
    assert ResourceOverrides.check('', known, ceilings).errors == []
}
Check.that('a valid block passes') {
    def r = check([PYWOMBAT: [memory: '64.GB', cpus: 4, time: '2.h', families: ['FAM*': [memory: '128.GB']]]])
    assert r == [errors: [], warnings: []]
}
Check.that('an unknown process, and a regex matching none, are errors') {
    assert check([PYWOMBATT: [cpus: 2]]).errors[0].contains('no process is called PYWOMBATT')
    assert check(['PY.*': [cpus: 2]]).errors == []
    assert check(['[': [cpus: 2]]).errors[0].contains('not a valid regex')
}
Check.that('an aliased process must be named as it runs') {
    def aliased = known.aliases.values().find()
    if (aliased) assert check([(aliased): [cpus: 2]]).errors[0].contains('runs as')
}
Check.that('bad values are errors: unitless memory, bad duration, fractional cpus, unknown field') {
    def errors = check([PYWOMBAT: [memory: 64, time: 'soon', cpus: 1.5, ram: '1.GB']]).errors
    assert errors.size() == 4
    assert errors.any { it.contains('.memory: expected a memory with its unit') }
    assert errors.any { it.contains('.time: expected a duration') }
    assert errors.any { it.contains('.cpus: expected a whole number') }
    assert errors.any { it.contains('.ram: unknown setting') }
}
Check.that('values above the ceilings are warnings, not errors') {
    def r = check([PYWOMBAT: [memory: '1.TB', cpus: 128, time: '72.h']])
    assert r.errors == [] && r.warnings.size() == 3
}
Check.that('a non-map block or entry is an error') {
    assert ResourceOverrides.check('PYWOMBAT', known, ceilings).errors[0].startsWith('resources must be a map')
    assert check([PYWOMBAT: '8.GB']).errors[0].contains('must be a map')
    assert check([PYWOMBAT: [families: ['F1': [families: [:]]]]]).errors[0].contains('belongs directly under a process')
}

Check.done()
