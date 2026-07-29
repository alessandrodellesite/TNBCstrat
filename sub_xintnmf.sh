#!/usr/bin/env bash
#SBATCH --job-name=xintnmf_sweep
#SBATCH --output=xintnmf_%j.log
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=64G

IMAGE_PATH="/mnt/pet/alex/XintNMF/xintnmf_image.sif"
DATA_DIR="/mnt/pet/alex/XintNMF/data"
OUT_BASE="/mnt/pet/alex/XintNMF/output"

MAX_PARALLEL=8      # 32 cpus / 4 cpus-per-run = 8 concurrent runs
export OMP_NUM_THREADS=4

run_one () {
  local K=$1
  local RUN=$2
  local OUT_DIR="${OUT_BASE}/rankselect_k${K}_run${RUN}"
  if [ -f "${OUT_DIR}/sample_factor.csv" ]; then
      echo "k=${K} run=${RUN} already done, skipping."
      return
  fi
  echo "Running k=${K} run=${RUN}..."
  singularity exec \
    -B /home/alexdull@prbb.org:/home/alexdull@prbb.org,/mnt/pet/alex:/mnt/pet/alex \
    "$IMAGE_PATH" \
    python /opt/X-intNMF/X-intNMF-run.py \
    --omics_input "${DATA_DIR}/rna_processed.tsv" "${DATA_DIR}/methylation_processed.tsv" "${DATA_DIR}/cnv_processed.tsv" \
    --output_dir "$OUT_DIR" \
    --output_format csv \
    --num_components "$K" \
    --graph_regularization 0 \
    --max_iter 1500 \
    --backend numpy \
    --gpu -1
}

job_count=0
for K in 2 3 4 5 6; do
  for RUN in $(seq 1 15); do
    run_one "$K" "$RUN" &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$MAX_PARALLEL" ]; then
      wait
      job_count=0
    fi
  done
done
wait
