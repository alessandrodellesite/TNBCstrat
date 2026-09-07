#!/usr/bin/env bash
#SBATCH --job-name=snf_multiomics
#SBATCH --output=snf_multiomics%j.log
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL  
#SBATCH --mail-user=alessandrodelle@vhio.net
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=1:00:00      
#SBATCH --mem=10G            

IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/r_image.sif"
SCRIPT_PATH="/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/Multiomics_SNF.R"

# singularity execution
srun singularity exec \
  -B /home/alessandrodelle@vhio.org:/home/alessandrodelle@vhio.org,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
  $IMAGE_PATH \
  Rscript $SCRIPT_PATH


