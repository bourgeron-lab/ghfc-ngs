# Cohort run state: `.ghfc-ngs.state.json`

Every cohort directory gets a hidden JSON file recording what the pipeline last did to it:

```
${data}/cohorts/<COHORT_NAME>/.ghfc-ngs.state.json
```

By convention that is the same directory that already holds `<COHORT_NAME>.params.yml`,
`<COHORT_NAME>.pedigree.tsv` and the cohort's merged outputs, so the record sits with the thing
it describes. It answers, without reading a Nextflow log: when did this cohort last run, did it
finish, which pedigree and which parameters produced these outputs, and how far along is each
step.

The file is written by the pipeline. The mechanics live in
[`lib/CohortState.groovy`](lib/CohortState.groovy); what gets recorded is assembled in
[`main.nf`](main.nf) under `COHORT RUN STATE`.

## Guarantees, and what is deliberately not guaranteed

Read these before trusting a field.

- **Written only by the pipeline.** Nothing else should edit it. It is not in `.gitignore`'s way
  because it lives in the data tree, not the repo.
- **Writes are atomic.** A new version is written to a temporary file in the same directory and
  moved into place with `ATOMIC_MOVE`, so a reader sees either the whole previous version or the
  whole new one — never a half-written file, even mid-run.
- **Concurrent runs are last-writer-wins.** Two runs against one cohort will not corrupt the
  file, but the second to finish overwrites the first's `last_run`. There is no locking.
- **Not written by `-stub-run` or `-preview`.** A stub run fabricates real files with fake
  content; a state file reporting those steps as complete would be worse than no file at all.
- **Not written when `cohort_name` is unset.** `params.cohort_name` has no default in
  `nextflow.config`, and without it there is no path. The run logs a warning and continues.
- **Not written for a run that fails before `cohort_name` and `data` are known** — in practice,
  only a missing `--data`.
- **A missing file means "not run since this feature shipped"**, not "never run". Absence is not
  evidence.
- **`status: "running"` is not a liveness signal.** There is no heartbeat. See
  [Status lifecycle](#status-lifecycle).
- **Percentages from a failed run can under-report.** See `outputs_may_be_incomplete`.

## Top-level schema

| Field | Type | Meaning |
| --- | --- | --- |
| `schema_version` | int | Currently `1`. **Check this before parsing.** |
| `cohort_name` | string | The cohort this file describes. |
| `last_run` | record | The most recent run. May still be `running`. |
| `last_successful_run` | record \| absent | The most recent run that finished successfully. Only ever replaced by another success, so a later failure never erases it. |
| `history` | array of records | Up to 10 **terminal** records, newest first. A `running` record is never in here; it is added once it reaches a terminal status. |

`last_run` is usually also `history[0]`. The exception is while a run is in flight, when
`last_run` is `running` and not yet in `history`.

## A record

Records are self-contained — they repeat `cohort_name` — so that a tool can collect records
across many cohorts into one list without losing track of which is which.

```json
{
  "run_id": "6f1a2b3c-4d5e-6f70-8192-a3b4c5d6e7f8",
  "status": "success",
  "started_at": "2026-09-21T09:52:03+02:00",
  "finished_at": "2026-09-21T18:43:59+02:00",
  "cohort_name": "EAGER",
  "pipeline": {
    "version": "1.0.0",
    "revision": "main",
    "commit_id": "8f94f13a2c9d4e6b1f0a7c3d5e2b8a4f6c1d9e0b",
    "repository": "https://github.com/bourgeron-lab/ghfc-ngs"
  },
  "pedigree": {
    "path": "/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/cohorts/EAGER/EAGER.pedigree.tsv",
    "sha256": "7a4ebfacdd6b9891e739e94b746799a3787cb12740ac46ccfac114859fa0b869",
    "families": 314,
    "individuals": 941
  },
  "params_file": {
    "path": "/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/cohorts/EAGER/EAGER.params.yml",
    "sha256": "fb6d11d8c0335d176533755eb7d926635d6d013aae7fb2463d2481613d4a248c"
  },
  "params_effective_sha256": "4f50615493b5f007beb00925c53982225215ffc7a6b5f214256ff3c83f7557b2",
  "steps_requested": [
    "alignment", "deepvariant_sample", "deepvariant_family",
    "annotation", "snvs_cohort", "wisecondorx", "wombat"
  ],
  "completion_measured": "after",
  "outputs_may_be_incomplete": false,
  "completion": {
    "alignment":          { "done": 900, "total": 941, "pct": 95.64 },
    "deepvariant_sample": { "done": 880, "total": 941, "pct": 93.52 },
    "deepvariant_family": { "done": 300, "total": 314, "pct": 95.54 },
    "annotation":         { "done": 298, "total": 314, "pct": 94.9 },
    "wombat":             { "done": 290, "total": 314, "pct": 92.36 },
    "wisecondorx":        { "done": 850, "total": 941, "pct": 90.33 },
    "ancestry":           null,
    "snvs_cohort":        { "done": 1, "total": 1, "pct": 100.0 },
    "extractor":          null
  }
}
```

### Fields

| Field | Type | Null when | Meaning |
| --- | --- | --- | --- |
| `run_id` | string | — | Nextflow's session UUID. Use it to match a `running` record to its terminal record, and as the key that de-duplicates `history`. It is also what `nextflow run -resume <run_id>` takes. |
| `status` | string | — | `running`, `success`, `failed` or `interrupted`. See below. |
| `started_at` | string | — | ISO-8601 with offset. Nextflow's `workflow.start`. |
| `finished_at` | string | **while `running`** | ISO-8601 with offset. Wall-clock duration is `finished_at - started_at`; it is not stored separately. |
| `cohort_name` | string | — | Repeated from the top level so the record stands alone. |
| `pipeline.version` | string | — | `manifest.version` from `nextflow.config`. |
| `pipeline.revision` | string | — | Branch or tag, e.g. `main`. |
| `pipeline.commit_id` | string | **local runs** | Full SHA of the commit that ran. Null when launched as `nextflow run main.nf` from a working copy rather than as a remote project. **This is the only field that identifies the code that actually ran**, which matters because `run_pipeline.sh` launches with `-latest`. |
| `pipeline.repository` | string | **local runs** | Git remote URL. |
| `pedigree.path` | string | — | The pedigree `main.nf` read (`params.pedigree`, else `${data}/cohorts/${cohort_name}/${cohort_name}.pedigree.tsv`). |
| `pedigree.sha256` | string | **file absent** | See [Checksums](#checksums). |
| `pedigree.families` | int | no pedigree parsed | Distinct FIDs. The denominator for family-level steps. |
| `pedigree.individuals` | int | no pedigree parsed | Distinct barcodes. The denominator for sample-level steps. |
| `params_file` | object | **no `-params-file` given** | `{path, sha256}` of the YAML. Recovered from the command line, since Nextflow does not expose it. |
| `params_effective_sha256` | string | rarely, if params will not serialise | SHA-256 of the resolved parameter map. See [Checksums](#checksums). |
| `steps_requested` | array | — | `params.steps` for this run. **Needed to interpret `completion`** — a step that was not requested was not measured. |
| `completion_measured` | string | — | `after` = the tree was re-scanned once the run finished. `before` = the numbers are from the start of the run (every `failed`-at-validation record, and every `running` record). |
| `outputs_may_be_incomplete` | bool | — | `true` on a failed run. Nextflow cancels in-flight `publishDir` copies when a run aborts, so outputs this run produced may not have landed before the counts were taken. |
| `completion` | object | **no plan was built** | Per-step progress. See below. |
| `reason` | string | absent unless set | Short human explanation, present on records written by a validation failure, e.g. `invalid steps requested: bogus`. |

### A failed record

A run stopped by a validation error carries the same fields, with these differences. Everything
not listed here (`run_id`, `started_at`, `finished_at`, `cohort_name`, `pipeline`, `pedigree`,
`params_file`, `params_effective_sha256`) is present and populated exactly as above.

```json
{
  "status": "failed",
  "steps_requested": ["alignment", "deepvariant_sample", "bogus"],
  "completion_measured": "before",
  "outputs_may_be_incomplete": false,
  "completion": null,
  "reason": "invalid steps requested: bogus"
}
```

Two things to handle when reading a failure:

- **`completion` is `null`** when the run stopped before a plan was built, which is the case for
  every validation error. A run that failed during execution *does* carry a `completion` block,
  with `outputs_may_be_incomplete` set.
- **`reason` is only present on validation failures.** A run that failed because a task errored
  has no `reason` — the Nextflow log is the place to look for that. Treat `reason` as absent
  rather than null, and fall back to `status` alone.


## Status lifecycle

```
                    validation passed, work starts
                                 │
                                 ▼
  (no file)  ──────────────► running ──────────┬──► success
                                 │             └──► failed
                                 │
                       run is killed and never
                       reaches the handler
                                 │
                                 ▼
                   stays "running" until the NEXT
                   run rewrites it to "interrupted"
                   and files it in history
```

| Status | Written when |
| --- | --- |
| `running` | Validation has passed and real work is starting. Not in `history`. |
| `success` | The run finished and Nextflow reported success. Also becomes `last_successful_run`. |
| `failed` | A validation error, or a task failure. Carries `reason` for the former. |
| `interrupted` | Never written by the run it describes. A later run found a `running` record from a *different* `run_id` and closed it out. |

**A `running` record does not mean a run is in progress.** Nextflow's `exit` is a `System.exit`
with no shutdown hook, so a SLURM walltime kill, a Ctrl-C, or a lost launch node leaves the
`running` record behind with nothing to close it. If `last_run.status` is `running` and
`started_at` is not recent, that run died. The next run against the cohort will move it to
`history` as `interrupted`.

## Reading `completion`

Each step is either `null` or `{done, total, pct}`, where `pct` is `done/total` as a percentage
rounded to two decimals.

**`null` means unmeasured, not 0%.** Two steps can be null:

- `ancestry` — when `ancestry` is not in `steps_requested`. The pipeline skips that whole part of
  the planning pass, so there is nothing to report.
- `extractor` — always. The pipeline has no on-disk completeness check for it.

Other things worth knowing before you quote a number:

- **`total` comes from the pedigree**, not from the pipeline's own "done + to do" counts. Those
  two are not a partition of the cohort: work is only scheduled when its prerequisites are met,
  so an individual blocked upstream appears in neither, and their sum understates the cohort.
- **`alignment` counts CRAMs.** A sample whose CRAM exists but whose coverage bedgraph is missing
  still counts as done — only the bedgraph is regenerated, never the alignment.
- **`snvs_cohort` is cohort-level**: `done` is 0 or 1 out of 1, not a count of families.
- **`completion_measured: "before"` means the numbers predate the run.** A validation failure and
  a `running` record both report the tree as it was at launch.

| Step | Denominator |
| --- | --- |
| `alignment`, `deepvariant_sample`, `wisecondorx` | `pedigree.individuals` |
| `deepvariant_family`, `annotation`, `wombat`, `ancestry` | `pedigree.families` |
| `snvs_cohort` | 1 |

## Checksums

Two different questions, two different fields.

**`pedigree.sha256` and `params_file.sha256`** — "has this file been edited since the run?"
Plain SHA-256 of the file's bytes. Reproduce them with:

```bash
shasum -a 256 /pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/cohorts/EAGER/EAGER.pedigree.tsv
```

`pedigree.sha256` covers the cohort pedigree that `main.nf` read. It does **not** cover the
per-family `{FID}.pedigree.tsv` files, which the pipeline regenerates from it.

**`params_effective_sha256`** — "was this run configured the same way?" SHA-256 of the resolved
parameter map, serialised as JSON with keys sorted. It catches what the file hash cannot:
`run_pipeline.sh` lets `--steps`, `--data`, `--pedigree`, `--ref` and `--ref-name` override the
YAML, and those change results without changing the file. It is opaque by design — it tells you
*that* something differs, never *what* — and it is **not comparable across pipeline versions**
that add, remove or rename a parameter, since a new default shifts the hash.

Neither checksum is verified on the next run. This file is a record, not a guard.

## Recipes

```bash
# Did the last run work, and when did this cohort last complete?
jq -r '.last_run.status' .ghfc-ngs.state.json
jq -r '.last_successful_run.finished_at // "never"' .ghfc-ngs.state.json
```

```bash
# Per-step progress as a table
jq -r '.last_run.completion | to_entries[]
       | "\(.key)\t\(if .value then "\(.value.done)/\(.value.total) (\(.value.pct)%)" else "not measured" end)"' \
   .ghfc-ngs.state.json | column -t -s $'\t'
```

```bash
# Fleet view: every cohort, its last status and its last success
find "$GHFC_NGS_COHORTS" -maxdepth 2 -name .ghfc-ngs.state.json -print0 \
  | xargs -0 jq -r '[.cohort_name, .last_run.status,
                     (.last_run.started_at // "-"),
                     (.last_successful_run.finished_at // "never")] | @tsv' \
  | sort | column -t -s $'\t'
```

```bash
# Has the pedigree been edited since the last successful run?
f=.ghfc-ngs.state.json
recorded=$(jq -r '.last_successful_run.pedigree.sha256' "$f")
current=$(shasum -a 256 "$(jq -r '.last_successful_run.pedigree.path' "$f")" | cut -d' ' -f1)
[ "$recorded" = "$current" ] && echo "pedigree unchanged" || echo "PEDIGREE CHANGED since last success"
```

```bash
# Runs that died without finishing (a stale 'running' is one of them)
jq -r 'select(.last_run.status == "running")
       | "\(.cohort_name) has an unfinished run started \(.last_run.started_at)"' .ghfc-ngs.state.json
```

## Extending the schema

New fields are additive: a consumer must ignore keys it does not know. Bump `schema_version`
only when an existing field changes meaning or disappears.

To add a field, extend `buildStateRecord` in [`main.nf`](main.nf). To change how records are
folded into the file — the history rules, the `interrupted` promotion, the atomic write — edit
[`lib/CohortState.groovy`](lib/CohortState.groovy), which has no Nextflow dependency and can be
exercised on its own.
