#!/usr/bin/env bash
#SBATCH --job-name=mofa2_boot
#SBATCH --output=logs/mofa2_boot_%A_%a.out
#SBATCH --array=1-1000%15
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL  
#SBATCH --mail-user=alessandrodelle@vhio.net
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2 
#SBATCH --time=45:00      
#SBATCH --mem=20G       

mkdir -p logs
IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/r_image.sif"
SCRIPT_PATH="/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/Bootstrap_mofa.R"
OUTDIR="/mnt/petasan_ccb/alessandro/SCANB/multiomics/bootstrap/mofa/mofa_bootstrap_k2"

srun singularity exec \
  --no-home \
  --env R_LIBS_USER="" \
  --env R_LIBS="" \
  --env HOME=/tmp \
  -B /home/alessandrodelle@vhio.org:/home/alessandrodelle@vhio.org,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
  $IMAGE_PATH \
  Rscript $SCRIPT_PATH "${SLURM_ARRAY_TASK_ID}" "${OUTDIR}"

