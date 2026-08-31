#!/usr/bin/env bash
#SBATCH --job-name=xintnmf_reg
#SBATCH --output=xintnmf_%j.log
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alessandrodelle@vhio.net
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=9:00:00
#SBATCH --mem=64G

## xintNMF with regularization 

IMAGE_PATH="/mnt/petasan_ccb/alessandro/SCANB/xint_image.sif"
DATA_DIR="/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/xintnmf_inputdata"
OUT_BASE="/mnt/petasan_ccb/alessandro/SCANB/multiomics/xintnmf_results"

MAX_PARALLEL=5      # one task per k, 4 cpus each = 20 cpus
export OMP_NUM_THREADS=4

run_one () {
  local K=$1
  local OUT_DIR="${OUT_BASE}/rankselect_k${K}_graphreg"
  if [ -f "${OUT_DIR}/sample_factor.csv" ]; then
      echo "k=${K} already done, skipping."
      return
  fi
  echo "Running k=${K}..."
  singularity exec \
    -B /home/alessandrodelle@vhio.org:/home/alessandrodelle@vhio.org,/mnt/petasan_ccb/alessandro:/mnt/petasan_ccb/alessandro \
    "$IMAGE_PATH" \
    python /opt/X-intNMF/X-intNMF-run.py \
    --omics_input "${DATA_DIR}/rna.tsv" "${DATA_DIR}/methylation.tsv" "${DATA_DIR}/cnv.tsv" \
    --interaction_input "${DATA_DIR}/interaction_rna_methylation.tsv" "${DATA_DIR}/interaction_rna_cnv.tsv" \
    --output_dir "$OUT_DIR" \
    --output_format csv \
    --num_components "$K" \
    --graph_regularization 1 \
    --max_iter 5000 \
    --backend numpy \
    --gpu -1
}

job_count=0
for K in 2 3 4 5 6; do
  run_one "$K" &
  job_count=$((job_count + 1))
  if [ "$job_count" -ge "$MAX_PARALLEL" ]; then
    wait
    job_count=0
  fi
done
wait
