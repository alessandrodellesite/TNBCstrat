#!/usr/bin/env bash
#SBATCH --job-name=snf4_boot
#SBATCH --output=logs/snf4_boot%A_%a.out
#SBATCH --array=1-1000%25
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL  
#SBATCH --mail-user=alessandrodelle@vhio.net
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2 
#SBATCH --time=15:00      
#SBATCH --mem=3G       

mkdir -p logs
#mkdir -p /mnt/petasan_ccb/alessandro/SCANB/basilisk_cache #for mofa only

IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/r_image.sif"
SCRIPT_PATH="/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/Bootstrap_snf.R"
OUTDIR="/mnt/petasan_ccb/alessandro/SCANB/multiomics/bootstrap/snf/snf_bootstrap_k4"


#for mofa:
#srun singularity exec \
#  --no-home \
#  --env XDG_CACHE_HOME=/mnt/petasan_ccb/alessandro/SCANB/basilisk_cache \
#  -B /home/alessandrodelle@vhio.org/ondemand/TNBCstrat:/home/alessandrodelle@vhio.org/ondemand/TNBCstrat,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
#  $IMAGE_PATH \
#  Rscript $SCRIPT_PATH "${SLURM_ARRAY_TASK_ID}" "${OUTDIR}"

srun singularity exec \
  -B /home/alessandrodelle@vhio.org:/home/alessandrodelle@vhio.org,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
  $IMAGE_PATH \
  Rscript $SCRIPT_PATH "${SLURM_ARRAY_TASK_ID}" "${OUTDIR}"
