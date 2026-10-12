#!/usr/bin/env bash

# Runs the self-contained *_test scripts that sit next to the scripts they cover, and reports
# which failed. With a directory argument, only the tests under it.
#
#   tests/run_script_tests.sh               # all of them
#   tests/run_script_tests.sh scripts       # the operator tools only, e.g. on the cluster's python

set -uo pipefail

cd "$(dirname "$0")/.." || exit 1

tests=()
while IFS= read -r t; do tests+=("${t}"); done < <(find "${@:-.}" -type f -name '*_test' -not -path '*/.git/*' -not -path '*/work/*' | sed 's#^\./##' | sort)
if (( ${#tests[@]} == 0 )); then
    echo "no *_test scripts found under ${*:-.}" >&2
    exit 1
fi

failed=()
for t in "${tests[@]}"; do
    echo "=== ${t}"
    if ! "./${t}"; then failed+=("${t}"); fi
done

echo
echo "${#tests[@]} test scripts, ${#failed[@]} failed"
for t in ${failed[@]+"${failed[@]}"}; do echo "  FAILED  ${t}"; done
(( ${#failed[@]} == 0 ))
