#!/bin/bash

module load graalvm
module load apptainer
module load graphviz
module load nextflow

set -euo pipefail

ulimit -v unlimited
ulimit -Sn 65536
ulimit -u 65536
export OPENBLAS_NUM_THREADS=1

# Function to display usage
usage() {
    cat << EOF
GHFC WGS Family-based Variant Calling Pipeline

USAGE:
    $0 [OPTIONS]

OPTIONS:
    --profile PROFILE           Nextflow profile(s) to use (default: slurm,apptainer)
    --config CONFIG             Additional Nextflow config file
    --params-file FILE          Parameters file (YAML format)
    --work-dir DIR              Nextflow work directory (default: work)
    --data DIR                  Data directory
    --scratch DIR               Scratch directory  
    --pedigree FILE             Pedigree file (TSV format)
    --ref FILE                  Reference genome file
    --ref-name NAME             Reference genome name
    --steps "step1,step2"       Pipeline steps to run (alignment,deepvariant,family_calling)
    --resume                    Resume previous run
    --dry-run                   Show what would be executed
    --stub-run                  Run in stub mode (for testing)
    -h, --help                  Show this help message

EXAMPLES:
    # Run full pipeline
    $0 --params-file params.yml

    # Run with specific steps
    $0 --params-file params.yml --steps "deepvariant,family_calling"

    # Resume previous run
    $0 --params-file params.yml --resume

EOF
}

# Default values
PROFILE="slurm,apptainer"
CONFIG=""
RESUME=""
WORK_DIR="work"
PARAMS_FILE=""
DATA=""
SCRATCH=""
PEDIGREE=""
REF=""
REF_NAME=""
STEPS=""
DRY_RUN=""
STUB_RUN=""
EXTRA_ARGS=""

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

# Check if Nextflow is available
if ! command -v nextflow &> /dev/null; then
    echo "ERROR: Nextflow is not available. Please load the nextflow module."
    exit 1
fi

# Construct the command
CMD="nextflow run main.nf"
CMD="$CMD -profile $PROFILE"
CMD="$CMD -work-dir $WORK_DIR"

# Add optional parameters
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
