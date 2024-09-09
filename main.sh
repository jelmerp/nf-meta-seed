#!/bin/bash
#SBATCH --account=PAS2693
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nf-meta-%j.out
set -euo pipefail

# Load the Nextflow Conda environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/nextflow

# Constants
WORKFLOW=/fs/ess/PAS2693/jelmer/meta_pipeline
WORKDIR=/fs/scratch/PAS2693/jelmer/nf-meta
OUTDIR=results/nf-meta

# Report
echo
date
echo -e "\n# Starting Nextflow run with Nextflow base call:"
echo "nextflow run $WORKFLOW -ansi-log false -resume -work-dir $WORKDIR" 
echo -e "\n# ... and with pipeline parameters:"
echo "--outdir $OUTDIR $*"
echo -e "\n==========================================\n"

# Run the workflow
nextflow run $WORKFLOW \
    -ansi-log false \
    -resume \
    --outdir "$OUTDIR" \
    -work-dir "$WORKDIR" \
    "$@"

# Report
echo
date
