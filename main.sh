#!/bin/bash
#SBATCH --account=PAS2693
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nf-meta-%j.out
set -euo pipefail

# Load the Nextflow Conda environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/nextflow

# Script options and args
reads=
outdir=
more_opts=()
while [ "$1" != "" ]; do
    case "$1" in
        -i | --reads )      shift; reads=$1 ;;
        -o | --outdir )     shift; outdir=$1 ;;
        * )                 more_opts+=("$1") ;;
    esac
    shift
done
[[ -z "$reads" ]] && echo "ERROR: Please use --reads <reads> to specify your FASTQ files" && exit 1
[[ -z "$outdir" ]] && echo "ERROR: Please use --outdir <dir> to specify your output dir" && exit 1

# Constants
WORKFLOW=/fs/ess/PAS2693/jelmer/workflows/nf-meta
WORKDIR=/fs/scratch/PAS2693/jelmer/nf-meta

# Report
echo
date
echo -e "\n# Starting Nextflow run with Nextflow base call:"
echo "nextflow run $WORKFLOW -ansi-log false -resume -work-dir $WORKDIR" 
echo -e "\n# ... and with pipeline parameters:"
echo "--outdir $outdir --reads $reads ${more_opts[*]}"
echo -e "\n==========================================\n"

# Run the workflow
nextflow run $WORKFLOW \
    -ansi-log false \
    -resume \
    -work-dir "$WORKDIR" \
    --reads "$reads" \
    --outdir "$outdir" \
    ${more_opts[*]}

# Report
echo
date
