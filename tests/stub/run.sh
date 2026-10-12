#!/usr/bin/env bash

# Runs the whole pipeline on a fake cohort (make_fixture.sh), as CI does:
#
#   1. validation - parameter mistakes the run must refuse before doing anything, each with
#      its message (-preview, so nothing is written even if a check stops working);
#   2. fresh      - all 12 steps with -stub-run on the untouched cohort;
#   3. rerun      - the same again, on the tree the fresh run published.
#
# A pass fails on any FAILED task, and when the tasks it ran - name and tag, sorted - differ from
# expected/<pass>.tasks.txt. The rerun is the planner's own test: anything it schedules is
# something it did not recognise as done. Only the extractor, which has no completeness check,
# is expected there.
#
# A change that means to alter which tasks run regenerates the lists; review their diff:
#
#   tests/stub/run.sh [--update]
#
# Uses `nextflow` from PATH (NEXTFLOW to override) and honours NXF_VER through the launcher.
# Works in STUB_OUT, a fresh temporary directory by default, kept for inspection.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "${here}/../.." && pwd)
nextflow=${NEXTFLOW:-nextflow}
update=0
[[ "${1:-}" == --update ]] && update=1

out=${STUB_OUT:-$(mktemp -d "${TMPDIR:-/tmp}/ghfc-ngs-stub.XXXXXX")}
rm -rf "${out}"
mkdir -p "${out}"
echo "stub runs in ${out}"

"${here}/make_fixture.sh" "${out}/fixture"
base="${out}/fixture/base.yml"
fixture="${out}/fixture"

fail=0

# A params file: the fixture's base plus the given YAML lines
params() {
    local name=$1; shift
    { cat "${base}"; printf '%s\n' "$@"; } > "${out}/${name}.yml"
    echo "${out}/${name}.yml"
}

valid_scenario=(
    'cohort_name: "STUB"'
    'steps: ["alignment", "deepvariant_sample", "deepvariant_family", "annotation", "wombat", "wombat_cohort", "snvs_cohort", "wisecondorx", "extractor", "ancestry", "roh", "expansionhunter"]'
    "ancestry_reference: \"${fixture}/refs/bundle\""
    'expansionhunter_catalog_name: "eh_stub"'
)

# --- 1. validation ------------------------------------------------------------------------

# name | expected message | params lines (replacing the valid scenario's)
expect_refusal() {
    local name=$1 message=$2; shift 2
    local file log status=0
    file=$(params "validation_${name}" "$@")
    log="${out}/validation_${name}.out"
    (cd "${out}" && "${nextflow}" run "${root}/main.nf" -preview -params-file "${file}" \
        -c "${here}/ci.config" > "${log}" 2>&1) || status=$?
    if (( status == 0 )); then
        echo "  FAIL  validation ${name}: the run was accepted"; fail=1
    elif ! grep -qF -- "${message}" "${log}"; then
        echo "  FAIL  validation ${name}: refused, but without '${message}' (see ${log})"; fail=1
    else
        echo "  ok    validation ${name}"
    fi
}

echo "=== validation"
v=("${valid_scenario[@]}")
expect_refusal no_cohort_name "cohort_name parameter is required" "${v[1]}" "${v[2]}" "${v[3]}"
expect_refusal invalid_step "Invalid steps specified: bogus" "${v[0]}" 'steps: ["alignment", "bogus"]' "${v[2]}" "${v[3]}"
expect_refusal roh_without_tier "needs the reference bundle's roh/ tier" "${v[0]}" 'steps: ["roh"]' \
    "ancestry_reference: \"${fixture}/refs\"" "${v[3]}"
expect_refusal eh_bad_catalog_name "may only hold letters" "${v[0]}" 'steps: ["expansionhunter"]' "${v[2]}" \
    'expansionhunter_catalog_name: "eh/stub"'
expect_refusal resources_typo "no process is called PYWOMBATT" "${v[@]}" 'resources: {PYWOMBATT: {cpus: 2}}'

# --- 2. and 3. stub runs ------------------------------------------------------------------

file=$(params stub "${valid_scenario[@]}")
mkdir -p "${out}/launch"

stub_pass() {
    local pass=$1 status=0
    local log="${out}/${pass}.out" trace="${out}/${pass}.trace.txt" tasks="${out}/${pass}.tasks.txt"
    echo "=== ${pass}"
    (cd "${out}/launch" && "${nextflow}" run "${root}/main.nf" -stub-run -params-file "${file}" \
        -c "${here}/ci.config" -with-trace "${trace}" > "${log}" 2>&1) || status=$?
    if (( status != 0 )); then
        echo "  FAIL  ${pass}: nextflow exited ${status}"; tail -n 30 "${log}" | sed 's/^/        /'
        fail=1; return
    fi
    # The trace's name column: PROCESS (tag), with the workflow prefix
    if [[ -f "${trace}" ]]; then
        awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
            { print $col["name"] "\t" $col["status"] }' "${trace}" | LC_ALL=C sort > "${out}/${pass}.status.txt"
    else
        : > "${out}/${pass}.status.txt"
    fi
    cut -f1 "${out}/${pass}.status.txt" > "${tasks}"
    if grep -qv $'\tCOMPLETED$' "${out}/${pass}.status.txt"; then
        echo "  FAIL  ${pass}: tasks did not complete:"
        grep -v $'\tCOMPLETED$' "${out}/${pass}.status.txt" | sed 's/^/        /'
        fail=1
    fi
    local expected="${here}/expected/${pass}.tasks.txt"
    if (( update )); then
        cp "${tasks}" "${expected}"
        echo "  updated ${expected#"${root}"/} ($(wc -l < "${tasks}" | tr -d ' ') tasks)"
    elif ! diff -u "${expected}" "${tasks}" > "${out}/${pass}.diff"; then
        echo "  FAIL  ${pass}: the tasks run differ from ${expected#"${root}"/}:"
        sed 's/^/        /' "${out}/${pass}.diff"
        echo "        If that is intended, run tests/stub/run.sh --update and commit the lists."
        fail=1
    else
        echo "  ok    ${pass}: $(wc -l < "${tasks}" | tr -d ' ') tasks, as expected"
    fi
}

stub_pass fresh
stub_pass rerun

if (( fail )); then
    echo "stub runs FAILED - logs in ${out}"
    exit 1
fi
echo "stub runs passed"
