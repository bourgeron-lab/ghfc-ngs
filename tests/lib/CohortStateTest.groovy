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

Check.done()
