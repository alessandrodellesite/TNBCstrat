#!/usr/bin/env bash
#SBATCH --job-name=single_nmf
#SBATCH --output=single_nmf%j.log
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL  
#SBATCH --mail-user=alessandrodelle@vhio.net
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=32:00:00      
#SBATCH --mem=16G            

IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/r_image_new.sif"
SCRIPT_PATH="/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/NMF_all_singleomics.R"

# singularity execution
srun singularity exec \
  -B /home/alessandrodelle@vhio.org:/home/alessandrodelle@vhio.org,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
  $IMAGE_PATH \
  Rscript $SCRIPT_PATH

