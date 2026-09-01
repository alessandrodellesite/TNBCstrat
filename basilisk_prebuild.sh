#!/usr/bin/env bash
#SBATCH --job-name=basilisk_prebuild
#SBATCH --output=logs/basilisk_prebuild.out
#SBATCH --partition=short
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

mkdir -p /mnt/petasan_ccb/alessandro/SCANB/basilisk_cache
IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/r_image.sif"

singularity exec \
  --no-home \
  --env XDG_CACHE_HOME=/mnt/petasan_ccb/alessandro/SCANB/basilisk_cache \
  -B /home/alessandrodelle@vhio.org/ondemand/TNBCstrat:/home/alessandrodelle@vhio.org/ondemand/TNBCstrat,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
  $IMAGE_PATH \
  Rscript -e '
    library(MOFA2)
    m <- create_mofa(list(a = matrix(rnorm(100), 10, 10), b = matrix(rnorm(100), 10, 10)))
    m <- prepare_mofa(m)
    run_mofa(m, outfile = tempfile(fileext = ".hdf5"), use_basilisk = TRUE)
  '
