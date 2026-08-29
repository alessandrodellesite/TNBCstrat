#!/usr/bin/env bash
#SBATCH --job-name=methyl_adjust
#SBATCH --output=methyl_adjust%j.log
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL  
#SBATCH --mail-user=alessandrodelle@vhio.net
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --time=24:00      
#SBATCH --mem=128G            

IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/r_image.sif"
SCRIPT_PATH="/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/methyl_preprocess_1.R"

# singularity execution
srun singularity exec \
  -B /home/alessandrodelle@vhio.org:/home/alessandrodelle@vhio.org,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
  $IMAGE_PATH \
  Rscript $SCRIPT_PATH


