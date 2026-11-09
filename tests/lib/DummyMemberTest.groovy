// DummyMember decides which pedigree rows are never samples: a miss sends a dummy to alignment,
// where the run stops on a sample with no input.

Check.that('a dummy is told by its prefix') {
    assert DummyMember.isDummy('DUMMY_B00IDL0_F')
    assert DummyMember.isDummy('DUMMY_B00INOK_MF')
}
Check.that('a barcode, an unknown parent and nothing are not dummies') {
    assert !DummyMember.isDummy('B00IDL0')
    assert !DummyMember.isDummy('0')
    assert !DummyMember.isDummy(null)
}
Check.that('the prefix is case-sensitive and must lead') {
    assert !DummyMember.isDummy('dummy_B00IDL0_F')
    assert !DummyMember.isDummy('X_DUMMY_B00IDL0_F')
}

Check.done()
