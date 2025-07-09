#!/bin/bash
#SBATCH --account=PAS2693
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nfmeta-%j.out

# Load the Nextflow Conda environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/conda/nextflow-25.04

# Constants
WORKFLOW=jelmerp/nf-meta-seed
export NXF_SINGULARITY_CACHEDIR=~/containers

# Defaults
resume=true && resume_opt="-resume"   # Use Nextflow '-resume' option
local_wf=false                        # Use workflow from Github Repo

# Script options and args
params_file=
more_opts=()
while [ "$1" != "" ]; do
    case "$1" in
        --params )          shift; params_file=$1 ;;
        --restart )         resume=false ;;
        --local_wf )        shift; local_wf=true;WORKFLOW=$1 ;;
        * )                 more_opts+=("$1") ;;
    esac
    shift
done

# Strict Bash settings
set -euo pipefail

# Check options
[[ -z "$params_file" ]] && echo "ERROR: Please use --params <params-file> to specify your parameter file" && exit 1
[[ ! -f "$params_file" ]] && echo "ERROR: Parameter file $params_file does not exist" && exit 1
[[ "$resume" == false ]] && resume_opt=

# Report
echo
date

# Pull the latest version of the workflow
if [[ "$local_wf" == false ]]; then
    echo "# Getting the latest version of the workflow:"
    nextflow pull $WORKFLOW
    echo
fi

# Report
echo "# Starting Nextflow nf-meta run with the following command:"
echo "nextflow run $WORKFLOW -params-file $params_file -ansi-log false $resume_opt ${more_opts[*]}"
echo -e "\n==========================================\n"

# Run the workflow
nextflow run $WORKFLOW \
    -params-file "$params_file" \
    -ansi-log false \
    "$resume_opt" \
    ${more_opts[*]}

# Report
echo
date
