#!/bin/bash

# Nextflow cannot fetch a pipeline through the cluster proxy ("Text must not be null or empty"),
# and the login node reaches GitHub without it. A compute node is the opposite: GitHub only
# through the proxy. So a --submit job keeps it - it runs from code pulled at submit time and
# never fetches any, and Apptainer may still need the proxy to pull an image.
if [[ -z "${GHFC_NGS_IN_JOB:-}" ]]; then
    unset HTTP_PROXY https_proxy http_proxy HTTPS_PROXY
fi

module load graalvm/ce-java23-23.0.1
module load apptainer
module load graphviz
# module load nextflow

set -euo pipefail

ulimit -v unlimited
ulimit -Sn 65536
ulimit -u 65536
export OPENBLAS_NUM_THREADS=1
export NXF_ASSETS="${NXF_ASSETS:-$HOME/.nextflow/assets}"

# Function to display usage
usage() {
    cat << EOF
GHFC WGS Family-based Variant Calling Pipeline

USAGE:
    $PROG_NAME COHORT [OPTIONS]
    $PROG_NAME --params-file FILE [OPTIONS]

ARGUMENTS:
    COHORT                      Cohort name, shorthand for
                                --params-file cohorts/COHORT/COHORT.params.yml.
                                Looked up under ./cohorts first, then under
                                \$GHFC_NGS_COHORTS
                                (default: $COHORTS_ROOT)
                                Must be the first argument. Nextflow is then launched
                                from \$GHFC_NGS_RUNS/COHORT, one directory per cohort
                                (default: runs/ beside the cohorts/ directory the cohort
                                was found in), so each cohort has its own .nextflow.log,
                                history, reports and resume target.

OPTIONS:
    --profile PROFILE           Nextflow profile(s) to use (default: slurm,apptainer)
    --config CONFIG             Additional Nextflow config file
    --params-file FILE          Parameters file (YAML format), overrides COHORT
    --work-dir DIR              Nextflow work directory (default: work)
    --data DIR                  Data directory
    --scratch DIR               Scratch directory  
    --pedigree FILE             Pedigree file (TSV format)
    --ref FILE                  Reference genome file
    --ref-name NAME             Reference genome name
    --steps "step1,step2"       Pipeline steps to run (alignment,deepvariant,family_calling)
    --clean-stale-families      Delete the outputs of families whose pedigree gained members
                                that were never called, so they are re-called in full. Only
                                touches families that can actually be rebuilt from the data
                                on disk and the steps requested; reports the rest. Combine
                                with --dry-run to see what it would remove.
    --migrate                   Run migration workflow instead of main pipeline
    --resume [SESSION]          Resume a previous run. With a COHORT and no SESSION, the
                                cohort's last run (run_id in its .ghfc-ngs.state.json);
                                otherwise the given session ID or run name, or the last
                                run of the launch directory.
    --submit                    Run Nextflow itself as a Slurm job on ghfc, named
                                ghfc-ngs.PROJECT.COHORT (PROJECT: the last directory
                                of data:), instead of in this shell. Prints the job
                                ID and returns. Head job size: \$GHFC_NGS_HEAD_CPUS
                                (default 2), \$GHFC_NGS_HEAD_MEM (16G), JVM heap
                                \$GHFC_NGS_HEAD_HEAP (12g). The pipeline is pulled here,
                                and the job runs that checkout.
    --here                      Launch from the current directory, not \$GHFC_NGS_RUNS/COHORT

ENVIRONMENT:
    GHFC_NGS_REVISION           Branch, tag or commit to run instead of main
    --dry-run                   Show what would be executed
    --stub-run                  Run in stub mode (for testing)
    -h, --help                  Show this help message

EXAMPLES:
    # Run full pipeline for a cohort
    $PROG_NAME CANDY_mpx

    # Run full pipeline with an explicit parameters file
    $PROG_NAME --params-file params.yml

    # Run migration workflow (one-time, for legacy files)
    $PROG_NAME --migrate --params-file params.yml

    # Run with specific steps
    $PROG_NAME --params-file params.yml --steps "deepvariant,family_calling"

    # Resume previous run
    $PROG_NAME --params-file params.yml --resume

    # Run a cohort as a Slurm job, then resume its last run the same way
    $PROG_NAME CANDY_mpx --submit
    $PROG_NAME CANDY_mpx --submit --resume

    # Stop a submitted run cleanly: Nextflow cancels its own tasks on SIGTERM
    scancel --signal=TERM --batch --name=ghfc-ngs.GHFC-GRCh38.CANDY_mpx

    # A cohort of another project
    GHFC_NGS_COHORTS=/pasteur/helix/projects/ghfc_wgs/WES/SPARK-GRCh38/cohorts $PROG_NAME test --submit

    # See what a stale-family clean would delete, then do it
    $PROG_NAME CANDY_mpx --clean-stale-families --dry-run
    $PROG_NAME CANDY_mpx --clean-stale-families

EOF
}

# Function to abort with a message on stderr, before nextflow is launched
die() {
    echo "ERROR: $*" >&2
    exit 1
}

# Function to make a path absolute, since the run may be launched from another directory
abspath() {
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *)  printf '%s\n' "$PWD/$1" ;;
    esac
}

# Function to read last_run.run_id from a cohort's state file; prints nothing if unknown
last_run_id() {
    local state="$1/.ghfc-ngs.state.json"
    [[ -f "$state" ]] || return 0
    python3 -c 'import json, sys
try:
    print(json.load(open(sys.argv[1]))["last_run"]["run_id"] or "")
except Exception:
    pass' "$state" 2>/dev/null || true
}

# Function to resolve a cohort name to its parameters file; sets PARAMS_FILE
resolve_cohort_params() {
    local name="$1"
    local local_dir="$PWD/cohorts/$name"
    local root_dir="$COHORTS_ROOT/$name"

    [[ -f "$name" ]] && die "'$name' is a file, not a cohort name. Did you mean: --params-file $name"
    if [[ ! "$name" =~ ^[A-Za-z0-9._-]+$ ]]; then
        die "invalid cohort name '$name'. Expected letters, digits, '.', '_' or '-' only. To use a parameters file by path: --params-file $name"
    fi
    if [[ -f "$local_dir/$name.params.yml" ]]; then
        PARAMS_FILE="-params-file $local_dir/$name.params.yml"
        COHORT_DIR="$local_dir"
        return 0
    fi
    if [[ -f "$root_dir/$name.params.yml" ]]; then
        PARAMS_FILE="-params-file $root_dir/$name.params.yml"
        COHORT_DIR="$root_dir"
        return 0
    fi
    if [[ -d "$local_dir" || -d "$root_dir" ]]; then
        die "cohort '$name': directory found but no parameters file in it. Expected one of:
         $local_dir/$name.params.yml
         $root_dir/$name.params.yml"
    fi
    die "cohort '$name': not found. No cohort directory at:
         $local_dir
         $root_dir
       Check the name, set GHFC_NGS_COHORTS, or pass --params-file FILE."
}

# Default values
PROFILE="slurm,apptainer"
CONFIG=""
RESUME=""
WORK_DIR=""
PARAMS_FILE=""
DATA=""
SCRATCH=""
PEDIGREE=""
REF=""
REF_NAME=""
STEPS=""
CLEAN_STALE=""
MIGRATE=""
DRY_RUN=""
STUB_RUN=""
EXTRA_ARGS=""
COHORT=""
COHORT_DIR=""
RESUME_ID=""
SUBMIT=""
HERE=""
PROG_NAME="ghfc-ngs"
COHORTS_ROOT="${GHFC_NGS_COHORTS:-/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/cohorts}"
COHORTS_ROOT="${COHORTS_ROOT%/}"
# Unset: runs/ next to the cohorts/ directory the cohort was found in, resolved below, so each
# project (GHFC-GRCh38, SPARK-GRCh38...) keeps its own launch directories
RUNS_ROOT="${GHFC_NGS_RUNS:-}"
RUNS_ROOT="${RUNS_ROOT%/}"
# Every run used to be launched from here, so the Nextflow caches of runs from before the move
# to per-cohort launch directories are found under it
LEGACY_LAUNCH_DIR="${GHFC_NGS_LEGACY_LAUNCH:-/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38}"
HEAD_CPUS="${GHFC_NGS_HEAD_CPUS:-2}"
HEAD_MEM="${GHFC_NGS_HEAD_MEM:-16G}"
HEAD_HEAP="${GHFC_NGS_HEAD_HEAP:-12g}"
# Branch, tag or commit to run instead of the default branch, e.g. to try a branch out
REVISION="${GHFC_NGS_REVISION:-}"

# Kept verbatim for --submit, which hands the same command line to the Slurm job
ORIG_ARGS=("$@")

# No arguments at all: there is nothing to run
if [[ $# -eq 0 ]]; then
    usage >&2
    die "no arguments. Give a cohort name (e.g. $PROG_NAME CANDY_mpx) or --params-file FILE."
fi

# A bare first argument is a cohort name, shorthand for
# --params-file cohorts/<NAME>/<NAME>.params.yml (resolved after parsing)
if [[ "$1" != -* ]]; then
    COHORT="$1"
    shift
fi

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --config)
            CONFIG="--config $(abspath "$2")"
            shift 2
            ;;
        --params-file)
            [[ $# -ge 2 ]] || die "--params-file requires a FILE argument."
            [[ -f "$2" ]] || die "parameters file not found: $2"
            PARAMS_FILE="-params-file $(abspath "$2")"
            shift 2
            ;;
        --work-dir)
            WORK_DIR="$(abspath "$2")"
            shift 2
            ;;
        --data)
            DATA="--data $(abspath "$2")"
            shift 2
            ;;
        --scratch)
            SCRATCH="--scratch $(abspath "$2")"
            shift 2
            ;;
        --pedigree)
            PEDIGREE="--pedigree $(abspath "$2")"
            shift 2
            ;;
        --ref)
            REF="--ref $(abspath "$2")"
            shift 2
            ;;
        --ref-name)
            REF_NAME="--ref-name $2"
            shift 2
            ;;
        --steps)
            STEPS="--steps $2"
            shift 2
            ;;
        --clean-stale-families|--clean_stale_families)
            # Nextflow maps a hyphenated --foo-bar to params.fooBar, never to params.foo_bar,
            # so the flag has to be handed over in the spelling nextflow.config declares
            CLEAN_STALE="--clean_stale_families"
            shift
            ;;
        --migrate)
            MIGRATE="migrate.nf"
            shift
            ;;
        --resume)
            RESUME="-resume"
            # An optional session: a UUID, or a Nextflow run name such as happy_goldberg
            if [[ $# -ge 2 && "$2" != -* ]]; then
                [[ "$2" =~ ^([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}|[a-z]+_[a-z]+)$ ]] \
                    || die "--resume: '$2' is neither a session ID nor a run name"
                RESUME_ID="$2"
                shift
            fi
            shift
            ;;
        --submit)
            SUBMIT=1
            shift
            ;;
        --here)
            HERE=1
            shift
            ;;
        --dry-run)
            DRY_RUN="-preview"
            shift
            ;;
        --stub-run)
            STUB_RUN="-stub-run"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            EXTRA_ARGS="$EXTRA_ARGS $1"
            shift
            ;;
    esac
done

# Resolve the cohort shorthand, unless --params-file already gave one explicitly
if [[ -n "$COHORT" ]]; then
    if [[ -n "$PARAMS_FILE" ]]; then
        echo "WARNING: --params-file given; ignoring cohort name '$COHORT'." >&2
    else
        resolve_cohort_params "$COHORT"
    fi
fi

# The label of this run: the cohort name, else the params file's cohort_name. It names the
# Slurm job, which is how a live run is found again
LABEL="$COHORT"
if [[ -z "$LABEL" && -n "$PARAMS_FILE" ]]; then
    LABEL="$(sed -n 's/^cohort_name:[[:space:]]*["'"'"']\{0,1\}\([A-Za-z0-9._-]*\).*/\1/p' "${PARAMS_FILE#-params-file }" | head -1)"
fi
# The project: the directory the parameters' data: points at, e.g. GHFC-GRCh38 or SPARK-GRCh38,
# since that is where the outputs go. Cohort names repeat across projects - both have a
# `test` - so the project is part of the job name and of the tag on every task.
DATA_DIR="${DATA#--data }"
if [[ -z "$DATA_DIR" && -n "$PARAMS_FILE" ]]; then
    DATA_DIR="$(sed -n 's/^data:[[:space:]]*["'"'"']\{0,1\}\([^"'"'"' ]*\).*/\1/p' "${PARAMS_FILE#-params-file }" | head -1)"
fi
PROJECT=""
if [[ -n "$DATA_DIR" ]]; then
    PROJECT="$(basename "${DATA_DIR%/}")"
    PROJECT="${PROJECT//[^A-Za-z0-9_-]/_}"
fi
JOB_NAME="ghfc-ngs.${PROJECT:+$PROJECT.}${LABEL:-run}"
JOB_TAG="ghfc-ngs:${PROJECT:+$PROJECT/}${LABEL:-run}"

if [[ -z "$RUNS_ROOT" ]]; then
    if [[ -n "$COHORT_DIR" ]]; then
        RUNS_ROOT="$(dirname "$(dirname "$COHORT_DIR")")/runs"
    else
        RUNS_ROOT="$(dirname "$COHORTS_ROOT")/runs"
    fi
fi

# One launch directory per cohort, unless --here. Paths given on the command line were made
# absolute above, so moving does not change what they point at.
if [[ -n "$COHORT" && -z "$HERE" ]]; then
    LAUNCH_DIR="$RUNS_ROOT/$COHORT"
else
    LAUNCH_DIR="$PWD"
fi

# A second real run of the same cohort would race the first one on every output it publishes.
# Inside the submitted job, the job itself is the one squeue finds, so it is left out.
if [[ -n "$LABEL" && -z "$DRY_RUN" ]] && command -v squeue &> /dev/null; then
    # ghfc-ngs.<COHORT> too: the name jobs had before the project was part of it
    LIVE="$(squeue -h -u "$USER" -n "$JOB_NAME,ghfc-ngs.${LABEL}" -o %i 2>/dev/null | grep -vx "${SLURM_JOB_ID:-none}" | head -1 || true)"
    [[ -z "$LIVE" ]] || die "$JOB_NAME is already queued or running as Slurm job $LIVE. Stop it first: scancel --signal=TERM --batch $LIVE"
fi

# Hand the very same command line, minus --submit, to a Slurm job, and return
if [[ -n "$SUBMIT" ]]; then
    [[ -n "$LABEL" ]] || die "--submit needs a cohort name, or a params file that sets cohort_name"
    command -v sbatch &> /dev/null || die "--submit: sbatch is not available on $(hostname)"
    mkdir -p "$LAUNCH_DIR/.submit"
    STAMP="$(date +%Y%m%d-%H%M%S)"
    # The job runs a copy of this very script, so it runs what was just invoked rather than
    # whatever main is when it starts. Under the ghfc-ngs wrapper the script arrives via bash -c.
    RUNNER="$LAUNCH_DIR/.submit/run_pipeline.$STAMP.sh"
    if [[ -n "${BASH_EXECUTION_STRING:-}" ]]; then
        printf '%s\n' "$BASH_EXECUTION_STRING" > "$RUNNER"
    else
        cp "${BASH_SOURCE[0]}" "$RUNNER"
    fi
    JOB_ARGS=()
    for a in "${ORIG_ARGS[@]}"; do
        [[ "$a" == "--submit" ]] || JOB_ARGS+=("$a")
    done
    # The job starts from the directory the user submitted from, so the arguments mean exactly
    # what they meant here - relative paths, ./cohorts - and moves to the launch dir itself
    # bash -l, because `module` only exists in a login shell. exec all the way down - this
    # script ends in `exec nextflow`, which execs java - so the batch step IS Nextflow, and
    # `scancel --signal=TERM --batch` reaches it rather than a shell that would die first.
    # The job's node cannot fetch the pipeline (see the top of this script), so it is fetched
    # here, into a clone of the cohort's own. It is never shared with another job, and never
    # moved under a live run - the check above saw to that - and it pins the code the job runs.
    command -v nextflow &> /dev/null || die "Nextflow is not available. Please load the nextflow module."
    JOB_ASSETS="$LAUNCH_DIR/.nextflow-assets"
    # Always an explicit revision: a pull without -r updates whatever branch the clone is on,
    # so one submit of a branch would otherwise leave every later one on that branch
    JOB_REVISION="${REVISION:-main}"
    echo "Pulling bourgeron-lab/ghfc-ngs ($JOB_REVISION) into $JOB_ASSETS"
    NXF_ASSETS="$JOB_ASSETS" nextflow -q pull bourgeron-lab/ghfc-ngs -r "$JOB_REVISION" \
        || die "could not pull bourgeron-lab/ghfc-ngs"
    JOB_SCRIPT="$LAUNCH_DIR/.submit/job.$STAMP.sh"
    {
        echo '#!/bin/bash'
        echo "export NXF_OPTS=\"-Xms1g -Xmx$HEAD_HEAP\""
        echo "export GHFC_NGS_IN_JOB=1"
        printf 'export NXF_ASSETS=%q\n' "$JOB_ASSETS"
        printf 'export GHFC_NGS_REVISION=%q\n' "$JOB_REVISION"
        printf 'cd %q\n' "$PWD"
        printf 'exec bash -l %q%s\n' "$RUNNER" "$(printf ' %q' "${JOB_ARGS[@]}")"
    } > "$JOB_SCRIPT"
    mkdir -p "$LAUNCH_DIR"
    JOB_ID="$(sbatch --parsable \
        -p ghfc --qos=ghfc -c "$HEAD_CPUS" --mem="$HEAD_MEM" \
        -J "$JOB_NAME" --comment="$JOB_TAG" \
        -D "$LAUNCH_DIR" -o "$LAUNCH_DIR/ghfc-ngs.%j.log" \
        "$JOB_SCRIPT")" || die "sbatch failed"
    echo "Submitted $JOB_NAME as Slurm job $JOB_ID"
    echo "  launch dir : $LAUNCH_DIR"
    echo "  console    : $LAUNCH_DIR/ghfc-ngs.$JOB_ID.log"
    echo "  stop       : scancel --signal=TERM --batch $JOB_ID"
    exit 0
fi

mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

# --resume with a cohort and no session: the cohort's own last run, which a bare -resume in a
# directory shared by several cohorts would not be
if [[ -n "$RESUME" && -z "$RESUME_ID" && -n "$COHORT_DIR" ]]; then
    RESUME_ID="$(last_run_id "$COHORT_DIR")"
    [[ -z "$RESUME_ID" ]] || echo "Resuming $COHORT's last run, $RESUME_ID"
fi
# A run launched before per-cohort launch directories keeps its task cache in the old shared
# one. Nextflow would open an empty cache here without a word, and re-run everything.
if [[ "$RESUME_ID" =~ ^[0-9a-f-]{36}$ && ! -d ".nextflow/cache/$RESUME_ID" \
      && -d "$LEGACY_LAUNCH_DIR/.nextflow/cache/$RESUME_ID" ]]; then
    echo "Copying the cache of session $RESUME_ID from $LEGACY_LAUNCH_DIR"
    mkdir -p .nextflow/cache
    cp -a "$LEGACY_LAUNCH_DIR/.nextflow/cache/$RESUME_ID" ".nextflow/cache/$RESUME_ID"
fi
[[ -z "$RESUME_ID" ]] || RESUME="-resume $RESUME_ID"

# Check if Nextflow is available
if ! command -v nextflow &> /dev/null; then
    die "Nextflow is not available. Please load the nextflow module."
fi


# Construct the command
if [[ -n "$MIGRATE" ]]; then
    CMD="nextflow run -latest bourgeron-lab/ghfc-ngs/$MIGRATE"
elif [[ -n "${GHFC_NGS_IN_JOB:-}" ]]; then
    # The checkout pulled at submit time: no -latest, which would fetch. -r always, since
    # Nextflow refuses to run a checkout off the default branch without it ("stuck on revision")
    CMD="nextflow run bourgeron-lab/ghfc-ngs -r ${REVISION:-main}"
else
    CMD="nextflow run -latest bourgeron-lab/ghfc-ngs${REVISION:+ -r $REVISION}"
fi
CMD="$CMD -profile $PROFILE"

# Add optional parameters
[[ -n "$WORK_DIR" ]] && CMD="$CMD -work-dir $WORK_DIR"
[[ -n "$CONFIG" ]] && CMD="$CMD $CONFIG"
[[ -n "$PARAMS_FILE" ]] && CMD="$CMD $PARAMS_FILE"
[[ -n "$DATA" ]] && CMD="$CMD $DATA"
[[ -n "$SCRATCH" ]] && CMD="$CMD $SCRATCH"
[[ -n "$PEDIGREE" ]] && CMD="$CMD $PEDIGREE"
[[ -n "$REF" ]] && CMD="$CMD $REF"
[[ -n "$REF_NAME" ]] && CMD="$CMD $REF_NAME"
[[ -n "$STEPS" ]] && CMD="$CMD $STEPS"
[[ -n "$CLEAN_STALE" ]] && CMD="$CMD $CLEAN_STALE"
[[ -n "$RESUME" ]] && CMD="$CMD $RESUME"
[[ -n "$DRY_RUN" ]] && CMD="$CMD $DRY_RUN"
[[ -n "$STUB_RUN" ]] && CMD="$CMD $STUB_RUN"
[[ -n "$EXTRA_ARGS" ]] && CMD="$CMD $EXTRA_ARGS"

echo "Launch dir: $PWD"
echo "Executing: $CMD"
echo ""

# Execute the command. exec, so that Nextflow takes this shell's place: a signal sent to the
# runner - scancel --signal=TERM --batch, for a submitted run - then reaches Nextflow itself.
eval exec $CMD
