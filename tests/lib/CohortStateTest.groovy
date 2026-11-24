// The state file is read by other people's tooling while runs are in flight; COHORT_STATE.md
// is the contract these checks hold the code to.

import java.nio.file.Files

def run(String id, String status, Map extra = [:]) { [run_id: id, status: status, cohort_name: 'C'] + extra }

Check.that('a first running record is last_run and not yet history') {
    def s = CohortState.merge([:], run('r1', 'running'))
    assert s.schema_version == CohortState.SCHEMA_VERSION
    assert s.last_run.run_id == 'r1' && s.history == [] && !s.containsKey('last_successful_run')
    assert s.keySet().toList() == ['schema_version', 'cohort_name', 'last_run', 'history']
}
Check.that('a success closes the run, once, into history and last_successful_run') {
    def s = CohortState.merge(CohortState.merge([:], run('r1', 'running')), run('r1', 'success'))
    assert s.last_successful_run.run_id == 'r1'
    assert s.history*.run_id == ['r1']
    s = CohortState.merge(s, run('r1', 'success'))
    assert s.history*.run_id == ['r1']
}
Check.that('a failure does not replace the last successful run') {
    def s = CohortState.merge(CohortState.merge([:], run('r1', 'success')), run('r2', 'failed'))
    assert s.last_run.run_id == 'r2' && s.last_successful_run.run_id == 'r1'
    assert s.history*.run_id == ['r2', 'r1']
}
Check.that('a running record left by another run becomes interrupted history') {
    def s = CohortState.merge(CohortState.merge([:], run('r1', 'running')), run('r2', 'running'))
    assert s.history.size() == 1 && s.history[0].run_id == 'r1' && s.history[0].status == 'interrupted'
}
Check.that('history is capped, newest first') {
    def s = [:]
    (1..5).each { s = CohortState.merge(s, run("r${it}", 'success'), 3) }
    assert s.history*.run_id == ['r5', 'r4', 'r3']
}
Check.that('progress: null when unmeasured, null pct for an empty total, two decimals otherwise') {
    assert CohortState.progress(null, 3) == null
    assert CohortState.progress(0, 0) == [done: 0, total: 0, pct: null]
    assert CohortState.progress(1, 3).pct == 33.33d
}
Check.that('checksums: null for a missing file, stable for a string') {
    assert CohortState.sha256('/no/such/file') == null && CohortState.sha256(null) == null
    assert CohortState.sha256OfString('') == 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855'
}
Check.that('write then read round-trips, atomically, and leaves no temp file') {
    def dir = Files.createTempDirectory('cohort-state')
    try {
        def s = CohortState.merge([:], run('r1', 'success', [steps: ['roh']]))
        CohortState.write(dir, s)
        assert CohortState.read(dir) == s
        assert dir.toFile().list().toList() == [CohortState.FILE_NAME]
    }
    finally { dir.toFile().deleteDir() }
}
Check.that('a missing or corrupt state file reads as empty') {
    def dir = Files.createTempDirectory('cohort-state')
    try {
        assert CohortState.read(dir) == [:] && CohortState.read(null) == [:]
        dir.resolve(CohortState.FILE_NAME).toFile().text = '{ not json'
        assert CohortState.read(dir) == [:]
    }
    finally { dir.toFile().deleteDir() }
}

// --- last_run_by_step ---------------------------------------------------------------------

ALL = ['alignment', 'wombat', 'sampleqc', 'genomeqc', 'roh', 'extractor']
def planned(String id, String status, List steps, Map extra = [:]) {
    def completion = ALL.collectEntries { s -> [(s): s in steps && s != 'extractor' ? CohortState.progress(5, 10) : null] }
    run(id, status, [steps_requested: steps, completion: completion, started_at: "t-${id}".toString(),
                     pedigree: [sha256: "ped-${id}".toString()],
                     params_file: [path: 'p.yml', sha256: "par-${id}".toString()]] + extra)
}

Check.that('a QC-only run after a full one leaves every other step its own entry') {
    def s = CohortState.merge([:], planned('r1', 'success', ['alignment', 'wombat', 'sampleqc', 'genomeqc']))
    s = CohortState.merge(s, planned('r2', 'success', ['sampleqc', 'genomeqc']))
    def by = s[CohortState.BY_STEP]
    assert by.keySet().toList() == ['alignment', 'genomeqc', 'sampleqc', 'wombat']
    assert by.alignment.last.run_id == 'r1' && by.wombat.last_success.run_id == 'r1'
    assert by.sampleqc.last.run_id == 'r2' && by.genomeqc.last_success.run_id == 'r2'
    assert by.wombat.last == [run_id: 'r1', status: 'success', started_at: 't-r1', pedigree_sha256: 'ped-r1',
                              params_sha256: 'par-r1', completion: [done: 5, total: 10, pct: 50.0d]]
    assert s.keySet().toList() == ['schema_version', 'cohort_name', 'last_run', 'last_successful_run',
                                   'last_run_by_step', 'history']
}
Check.that('a failed run is a step\'s last run, never its last success') {
    def s = CohortState.merge([:], planned('r1', 'success', ['roh']))
    s = CohortState.merge(s, planned('r2', 'failed', ['roh'], [outputs_may_be_incomplete: true]))
    def roh = s[CohortState.BY_STEP].roh
    assert roh.last.run_id == 'r2' && roh.last.outputs_may_be_incomplete == true
    assert roh.last_success.run_id == 'r1'
}
Check.that('a running record, an interruption and a validation failure record nothing') {
    def s = CohortState.merge([:], planned('r1', 'success', ['roh']))
    def before = s[CohortState.BY_STEP]
    s = CohortState.merge(s, planned('r2', 'running', ['roh']))
    s = CohortState.merge(s, planned('r3', 'running', ['roh']))   // r2 is filed as interrupted
    s = CohortState.merge(s, run('r4', 'failed', [steps_requested: ['roh', 'bogus'], completion: null]))
    assert s.history*.status.contains('interrupted')
    assert s[CohortState.BY_STEP] == before && before.roh.last.run_id == 'r1'
}
Check.that('only requested steps with a completion block count; extractor, never measured, is kept') {
    def by = CohortState.merge([:], planned('r1', 'success', ['roh', 'extractor', 'bogus']))[CohortState.BY_STEP]
    assert by.keySet().toList() == ['extractor', 'roh']
    assert by.extractor.last.completion == null && by.extractor.last.run_id == 'r1'
}
Check.that('a file written before the per-step record is seeded from its history and last success') {
    def r0 = planned('r0', 'success', ['alignment', 'wombat'])
    def r1 = planned('r1', 'failed', ['alignment'])
    def r2 = planned('r2', 'success', ['sampleqc'])
    def old = [schema_version: 1, cohort_name: 'C', last_run: r2, last_successful_run: r2, history: [r2, r1]]
    def by = CohortState.merge(old, planned('r3', 'running', ['roh']))[CohortState.BY_STEP]
    assert by.alignment.last.run_id == 'r1' && !by.alignment.containsKey('last_success')
    assert by.sampleqc.last_success.run_id == 'r2' && !by.containsKey('roh')
    // A last success older than the whole history is folded in first, history over it
    old = [schema_version: 1, cohort_name: 'C', last_run: r1, last_successful_run: r0, history: [r1]]
    by = CohortState.merge(old, planned('r3', 'running', ['roh']))[CohortState.BY_STEP]
    assert by.wombat.last_success.run_id == 'r0' && by.wombat.last.run_id == 'r0'
    assert by.alignment.last.run_id == 'r1' && by.alignment.last_success.run_id == 'r0'
}
Check.that('a record written twice does not change the per-step record') {
    def s = CohortState.merge([:], planned('r1', 'success', ['roh']))
    assert CohortState.merge(s, planned('r1', 'success', ['roh'])) == s
}
Check.that('the per-step record round-trips through the file') {
    def dir = Files.createTempDirectory('cohort-state')
    try {
        def s = CohortState.merge(CohortState.merge([:], planned('r1', 'success', ['roh', 'wombat'])),
                                  planned('r2', 'failed', ['roh']))
        CohortState.write(dir, s)
        def back = CohortState.read(dir)
        assert back == s
        assert CohortState.merge(back, planned('r3', 'success', ['wombat']))[CohortState.BY_STEP].roh.last.run_id == 'r2'
    }
    finally { dir.toFile().deleteDir() }
}

Check.done()
