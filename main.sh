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
reads=
outdir=
more_opts=()
while [ "$1" != "" ]; do
    case "$1" in
        -i | --reads )      shift; reads=$1 ;;
        -o | --outdir )     shift; outdir=$1 ;;
        -w | -work-dir )    shift; workdir=$1 ;;
        -restart )          resume=false ;;
        * )                 more_opts+=("$1") ;;
    esac
    shift
done

# Strict Bash settings
set -euo pipefail

# Check options
[[ -z "$reads" ]] && echo "ERROR: Please use -i/--reads <reads> to specify your FASTQ files" && exit 1
[[ -z "$outdir" ]] && echo "ERROR: Please use -o/--outdir <dir> to specify your output dir" && exit 1
[[ "$resume" == false ]] && resume_opt=

# Report
echo
date
echo "# Starting Nextflow nf-meta run with the following command:"
echo "nextflow run $WORKFLOW --reads $reads --outdir $outdir -work-dir $workdir -ansi-log false -resume ${more_opts[*]}"
echo -e "\n==========================================\n"

# Create the output dir
mkdir -p "$outdir"/logs

# Run the workflow
nextflow run $WORKFLOW \
    --reads "$reads" \
    --outdir "$outdir" \
    -work-dir "$workdir" \
    -ansi-log false \
    $resume_opt \
    ${more_opts[*]}

# Report
echo
date
