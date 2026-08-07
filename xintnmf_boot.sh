#!/usr/bin/env bash
#SBATCH --output=logs/xintnmf_boot_%A_%a.out
#SBATCH --array=1-1000%20
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=30:00
#SBATCH --mem=20G

# K is passed in at submit time via --export, e.g.:
#   sbatch --job-name=xintnmf_boot_k2 --export=ALL,K=2 xintnmf_boot.sh
#   sbatch --job-name=xintnmf_boot_k3 --export=ALL,K=3 xintnmf_boot.sh
#   sbatch --job-name=xintnmf_boot_k4 --export=ALL,K=4 xintnmf_boot.sh
if [ -z "$K" ]; then
  echo "ERROR: K not set. Submit with: sbatch --export=ALL,K=<2|3|4> xintnmf_boot.sbatch"
  exit 1
fi

mkdir -p logs

IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/xint_image.sif"
SCRIPT_PATH="/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/xintnmf_bootstrap.py"
OUTDIR="/mnt/petasan_ccb/alessandro/SCANB/xintnmf_bootstrap_k${K}"

export OMP_NUM_THREADS=4


srun singularity exec \
  -B /home/alessandrodelle@vhio.org:/home/alessandrodelle@vhio.org,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
  $IMAGE_PATH \
  Rscript $SCRIPT_PATH "${SLURM_ARRAY_TASK_ID}" "${K}" "${OUTDIR}"
