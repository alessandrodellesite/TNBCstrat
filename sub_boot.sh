#!/usr/bin/env bash
#SBATCH --job-name=mofa_boot
#SBATCH --output=logs/mofa_boot_k3_%A_%a.out
#SBATCH --array=1-1000%20
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL  
#SBATCH --mail-user=alessandrodelle@vhio.net
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2 
#SBATCH --time=1:00:00      
#SBATCH --mem=20G       

mkdir -p logs 

IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/r_image.sif"
SCRIPT_PATH="/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/mofa_bootstrap.R"
OUTDIR="/mnt/petasan_ccb/alessandro/SCANB/mofa_bootstrap_k3"

# singularity execution - arguments after Rscript $SCRIPT_PATH are passed through to commandArgs() inside icluster_bootstrap.R
srun singularity exec \
  -B /home/alessandrodelle@vhio.org:/home/alessandrodelle@vhio.org,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
  $IMAGE_PATH \
  Rscript $SCRIPT_PATH "${SLURM_ARRAY_TASK_ID}" "${OUTDIR}"

