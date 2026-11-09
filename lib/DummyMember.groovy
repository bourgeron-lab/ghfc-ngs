/**
 * Dummy pedigree members: individuals who were never sequenced, written into the pedigree only to
 * link real members - the shared parents of MZ twins, the grandparents that join an aunt to a
 * parent. A dummy's ID starts with DUMMY_, then the barcode of a sequenced child it hangs from and
 * the path up to it: F father, M mother, MF the mother's father, FM the father's mother...
 * (DUMMY_B00IDL0_F, DUMMY_B00INOK_MF). A dummy has a row of its own, with its sex and its parents,
 * which may be dummies too; one named only as a child's father or mother is a founder.
 *
 * A dummy is never a sample. The planner leaves it out of the individuals and of its family's
 * members, so no step looks for its reads, calls or coverage, and it is not reported as a sample
 * without input. It stays in the pedigree file, so the expected kinships go through it: twins
 * with the same dummy parents are siblings.
 *
 * The scripts under modules/ that read pedigree rows (cnv_coverage, eh_family_table) apply the
 * same prefix; keep them in step with this class.
 */
class DummyMember {

    static final String PREFIX = 'DUMMY_'

    static boolean isDummy(String id) {
        return id != null && id.startsWith(PREFIX)
    }
}
