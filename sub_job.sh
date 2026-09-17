#!/usr/bin/env bash
#SBATCH --job-name=mofa_1
#SBATCH --output=mofa%j.log
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL  
#SBATCH --mail-user=alessandrodelle@vhio.net
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00      
#SBATCH --mem=6G            

IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/r_image.sif"
SCRIPT_PATH="/home/alessandrodelle@vhio.org/ondemand/TNBCstrat/Multiomics_mofa.R"

# singularity execution
srun singularity exec \
  -B /home/alessandrodelle@vhio.org:/home/alessandrodelle@vhio.org,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro,/mnt/petasan_ccb/juanra:/mnt/petasan_ccb/juanra \
  $IMAGE_PATH \
  Rscript $SCRIPT_PATH
