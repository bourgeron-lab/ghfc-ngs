#!/usr/bin/env bash

# Checks on what is tracked, rather than on what it says:
#   - every tracked file but a .nf that starts with a shebang is executable - the tests and
#     several processes run scripts by their path, and a lost exec bit only shows up on the
#     cluster;
#   - nothing tracked is over 5 MB, or is a run artefact (work/, .nextflow*, reports/,
#     the cohort state file, .DS_Store);
#   - no leftover merge-conflict markers.
#
#   tests/checks/repo_hygiene.sh

set -euo pipefail

cd "$(dirname "$0")/../.."

fail=0
bad() { printf '  - %s\n' "$1"; fail=1; }

echo "repo_hygiene:"

while IFS= read -r -d '' entry; do
    mode=${entry%% *}
    path=${entry#*$'\t'}
    [[ -f "${path}" && "${path}" != *.nf ]] || continue
    if [[ "$(head -c 2 "${path}")" == '#!' && "${mode}" != 100755 ]]; then
        bad "${path} starts with a shebang but is not executable in git (git update-index --chmod=+x)"
    fi
done < <(git ls-files -s -z)

while IFS= read -r -d '' path; do
    [[ -f "${path}" ]] || continue
    size=$(wc -c < "${path}" | tr -d ' ')
    if (( size > 5 * 1024 * 1024 )); then
        bad "${path} is $(( size / 1024 / 1024 )) MB: keep large files out of the repository"
    fi
done < <(git ls-files -z)

while IFS= read -r path; do
    bad "${path} is a run artefact and should not be tracked"
done < <(git ls-files | grep -E '(^|/)(work|\.nextflow|reports|results)/|(^|/)\.nextflow\.log|\.ghfc-ngs\.state|(^|/)\.DS_Store$' || true)

while IFS= read -r hit; do
    bad "merge-conflict marker at ${hit}"
done < <(git grep -nE '^(<<<<<<< |>>>>>>> |=======$)' -- ':!*.md' || true)

if (( fail )); then exit 1; fi
echo "  ok"
