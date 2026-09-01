#!/usr/bin/env bash
#SBATCH --job-name=mofa4_boot
#SBATCH --output=logs/mofa4_boot_%A_%a.out
#SBATCH --array=1-1000%15
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL  
#SBATCH --mail-user=alessandrodelle@vhio.net
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2 
#SBATCH --time=45:00      
#SBATCH --mem=10G       

mkdir -p logs
mkdir -p /mnt/petasan_ccb/alessandro/SCANB/basilisk_cache

IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/r_image.sif"
  SCRIPT_PATH="/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/Bootstrap_mofa_4.R"
OUTDIR="/mnt/petasan_ccb/alessandro/SCANB/multiomics/bootstrap/mofa/mofa_bootstrap_k4"

srun singularity exec \
  --no-home \
  --env XDG_CACHE_HOME=/mnt/petasan_ccb/alessandro/SCANB/basilisk_cache \
  -B /home/alessandrodelle@vhio.org/ondemand/TNBCstrat:/home/alessandrodelle@vhio.org/ondemand/TNBCstrat,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
  $IMAGE_PATH \
  Rscript $SCRIPT_PATH "${SLURM_ARRAY_TASK_ID}" "${OUTDIR}"
