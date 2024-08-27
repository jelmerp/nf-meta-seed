#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nf-meta-%j.out
set -euo pipefail

# Load the Nextflow Conda environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/nextflow

# Point to the main workflow file/dir
# Eventually, this should point to a GitHub repo
WORKFLOW=/fs/ess/PAS2693/jelmer/meta_pipeline

# Report
date
echo "# Starting Nextflow run with command:"
echo nextflow run $WORKFLOW -ansi-log false -resume "$@"
echo -e "==========================================\n"

# Run the workflow
nextflow run $WORKFLOW -ansi-log false -resume "$@"

# Report
echo
date
