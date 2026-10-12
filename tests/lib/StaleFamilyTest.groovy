// StaleFamily is the only code that decides what --clean_stale_families deletes.

def all_steps = ['alignment', 'deepvariant_sample', 'deepvariant_family']

Check.that('family outputs are the joint call, annotation and wombat files, under the family dir') {
    def paths = StaleFamily.familyOutputs('/d', 'FAM1', 'vep115', ['dnm.yaml', 'rare.yml'], false)
    assert paths.every { it.startsWith('/d/families/1/M/FAM1/') }
    assert '/d/families/1/M/FAM1/vcfs/FAM1.norm.bcf' in paths
    assert '/d/families/1/M/FAM1/FAM1.pedigree.tsv' in paths
    assert '/d/families/1/M/FAM1/vcfs/FAM1.rare.vep115.annotated.bcf' in paths
    assert '/d/families/1/M/FAM1/wombat/FAM1.rare.vep115.annotated.dnm.tsv' in paths
    assert '/d/families/1/M/FAM1/wombat/FAM1.rare.vep115.annotated.rare.tsv' in paths
    assert paths.size() == paths.toSet().size()
}
Check.that('family outputs never reach the groups that heal themselves') {
    def paths = StaleFamily.familyOutputs('/d', 'FAM1', 'vep115', ['rare.yaml'], true)
    assert !paths.any { it.contains('/ancestry/') || it.contains('/roh/') || it.contains('/str/') || it.contains('/extractor/') || it.contains('.dnm.') }
}
Check.that('wisecondorx files only when the step is in play') {
    assert !StaleFamily.familyOutputs('/d', 'FAM1', 'v', [], false).any { it.contains('wisecondorx') }
    assert StaleFamily.familyOutputs('/d', 'FAM1', 'v', [], true).count { it.contains('wisecondorx') } == 2
}
Check.that('cohort outputs are grouped by the steps that rebuild them') {
    def groups = StaleFamily.cohortOutputs('/d', 'C', 'v', ['dnm.yaml'], true)
    assert groups*.steps == [['annotation', 'snvs_cohort'], ['annotation', 'wombat', 'wombat_cohort']]
    assert groups[0].paths.contains('/d/cohorts/C/svs/wisecondorx/C_aberrations.bed')
    assert groups[1].paths == ['/d/cohorts/C/wombat/C.rare.v.annotated.dnm.results.tsv']
    assert StaleFamily.cohortOutputs('/d', 'C', 'v', [], false).size() == 1
}
Check.that('a member with a CRAM only needs calling') {
    def r = StaleFamily.recovery(['KID'], ['KID'] as Set, [] as Set, [:], all_steps - 'alignment')
    assert r.recoverable && r.needs_alignment == []
}
Check.that('a member with resolvable input needs alignment, and the step') {
    assert StaleFamily.recovery(['KID'], [] as Set, ['KID'] as Set, [:], all_steps).recoverable
    def r = StaleFamily.recovery(['KID'], [] as Set, ['KID'] as Set, [:], all_steps - 'alignment')
    assert !r.recoverable && r.needs_alignment == ['KID']
    assert r.reasons.any { it.contains('alignment step was not requested') }
}
Check.that('a member with no input data is never recoverable') {
    def r = StaleFamily.recovery(['KID'], [] as Set, [] as Set, [:], all_steps)
    assert !r.recoverable && r.reasons == ['no input data on disk for KID']
}
Check.that('an unindexed input CRAM is named as such') {
    def r = StaleFamily.recovery(['KID'], [] as Set, [] as Set, [KID: [path: '/in/KID.cram']], all_steps)
    assert !r.recoverable && r.reasons[0].contains("unindexed for KID: /in/KID.cram")
}
Check.that('a family cannot be re-called without both deepvariant steps') {
    ['deepvariant_family', 'deepvariant_sample'].each { step ->
        def r = StaleFamily.recovery(['KID'], ['KID'] as Set, [] as Set, [:], all_steps - step)
        assert !r.recoverable && r.reasons.any { it.contains(step) }
    }
    assert !StaleFamily.recovery([], [] as Set, [] as Set, [:], null).recoverable
}
Check.that('config names drop their YAML suffix and empties') {
    assert StaleFamily.configNames(['a.yaml', 'b.yml', 'c', null, '']) == ['a', 'b', 'c']
    assert StaleFamily.configNames(null) == []
}

Check.done()
