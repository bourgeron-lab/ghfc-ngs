#!/bin/bash

unset HTTP_PROXY https_proxy http_proxy HTTPS_PROXY

module load graalvm/ce-java23-23.0.1
module load apptainer
module load graphviz
# module load nextflow

set -euo pipefail

ulimit -v unlimited
ulimit -Sn 65536
ulimit -u 65536
export OPENBLAS_NUM_THREADS=1
export NXF_ASSETS=$HOME/.nextflow/assets

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
                                Must be the first argument.

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
    --migrate                   Run migration workflow instead of main pipeline
    --resume                    Resume previous run
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

EOF
}

# Function to abort with a message on stderr, before nextflow is launched
die() {
    echo "ERROR: $*" >&2
    exit 1
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
        return 0
    fi
    if [[ -f "$root_dir/$name.params.yml" ]]; then
        PARAMS_FILE="-params-file $root_dir/$name.params.yml"
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
MIGRATE=""
DRY_RUN=""
STUB_RUN=""
EXTRA_ARGS=""
COHORT=""
PROG_NAME="ghfc-ngs"
COHORTS_ROOT="${GHFC_NGS_COHORTS:-/pasteur/helix/projects/ghfc_wgs/WGS/GHFC-GRCh38/cohorts}"
COHORTS_ROOT="${COHORTS_ROOT%/}"

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
            CONFIG="--config $2"
            shift 2
            ;;
        --params-file)
            [[ $# -ge 2 ]] || die "--params-file requires a FILE argument."
            [[ -f "$2" ]] || die "parameters file not found: $2"
            PARAMS_FILE="-params-file $2"
            shift 2
            ;;
        --work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        --data)
            DATA="--data $2"
            shift 2
            ;;
        --scratch)
            SCRATCH="--scratch $2"
            shift 2
            ;;
        --pedigree)
            PEDIGREE="--pedigree $2"
            shift 2
            ;;
        --ref)
            REF="--ref $2"
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
        --migrate)
            MIGRATE="migrate.nf"
            shift
            ;;
        --resume)
            RESUME="-resume"
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

# Check if Nextflow is available
if ! command -v nextflow &> /dev/null; then
    die "Nextflow is not available. Please load the nextflow module."
fi


# Construct the command
if [[ -n "$MIGRATE" ]]; then
    CMD="nextflow run -latest bourgeron-lab/ghfc-ngs/$MIGRATE"
else
    CMD="nextflow run -latest bourgeron-lab/ghfc-ngs"
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
[[ -n "$RESUME" ]] && CMD="$CMD $RESUME"
[[ -n "$DRY_RUN" ]] && CMD="$CMD $DRY_RUN"
[[ -n "$STUB_RUN" ]] && CMD="$CMD $STUB_RUN"
[[ -n "$EXTRA_ARGS" ]] && CMD="$CMD $EXTRA_ARGS"

echo "Executing: $CMD"
echo ""

# Execute the command
eval $CMD
