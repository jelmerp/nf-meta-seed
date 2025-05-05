#!/bin/bash
#SBATCH --account=PAS2693
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nf-meta-%j.out

# Load the Nextflow Conda environment
module load miniconda3
conda activate /fs/ess/PAS0471/jelmer/conda/nextflow

# Constants
WORKFLOW=/fs/ess/PAS2693/jelmer/nf-meta

# Defaults
workdir=work
resume=true && resume_opt="-resume"

# Script options and args
params=
more_opts=()
while [ "$1" != "" ]; do
    case "$1" in
        --params )          shift; params=$1 ;;
        --workdir )         shift; workdir=$1 ;;
        --restart )         resume=false ;;
        * )                 more_opts+=("$1") ;;
    esac
    shift
done

# Strict Bash settings
set -euo pipefail

# Check options
[[ -z "$params" ]] && echo "ERROR: Please use --params <params-file> to specify your parameter file" && exit 1
[[ "$resume" == false ]] && resume_opt=

# Report
echo
date
echo "# Starting Nextflow nf-meta run with the following command:"
echo "nextflow run $WORKFLOW -params-file $params -work-dir $workdir -ansi-log false $resume_opt ${more_opts[*]}"
echo -e "\n==========================================\n"

# Run the workflow
nextflow run $WORKFLOW \
    -params-file "$params" \
    -work-dir "$workdir" \
    -ansi-log false \
    $resume_opt \
    ${more_opts[*]}

# Report
echo
date
