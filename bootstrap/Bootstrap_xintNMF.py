#!/usr/bin/env python3
"""
Fits a single X-intNMF model (fixed K, graph-regularized) on a subsample
(without replacement) of the full cohort. Intended to be launched once per
SLURM array task. Subsample composition is read from the shared,
pre-generated subsample list (boot_subsamples.json, corresponding 
boot_subsamples.rds) so that replicate #i uses the exact same set of patients 
across every benchmarked method.

Usage (called automatically by the sbatch script):
    python xintnmf_bootstrap.py <iter_id> <k> <outdir>
"""
import sys
import os
import json
import shutil
import subprocess
import tempfile
import datetime

import pandas as pd


def log(msg):
    print(f"[{datetime.datetime.now().isoformat(timespec='seconds')}] {msg}", flush=True)


def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: python xintnmf_bootstrap.py <iter_id> <k> <outdir>")

    iter_id = int(sys.argv[1])
    k = int(sys.argv[2])
    outdir = sys.argv[3]

    data_dir = "/mnt/petasan_ccb/alessandro/SCANB/multiomics/input_data/xintnmf_inputdata"
    shared_dir = "/mnt/petasan_ccb/alessandro/SCANB/bootstrap_shared"

    os.makedirs(outdir, exist_ok=True)
    final_out = os.path.join(outdir, f"iter_{iter_id}")
    done_marker = os.path.join(final_out, "sample_factor.csv")

    if os.path.exists(done_marker):
        log(f"iter={iter_id} k={k} already done, skipping.")
        return

    log(f"Iteration: {iter_id} | k={k} | Outdir: {final_out}")

    # Load shared subsample 
    subsample_path = os.path.join(shared_dir, "boot_subsamples.json")
    with open(subsample_path) as f:
        subsample_list = json.load(f)

    key = str(iter_id)
    if key not in subsample_list:
        sys.exit(
            f"iter_id {iter_id} not found in {subsample_path} "
            f"(max = {len(subsample_list)})"
        )

    sample_ids_sub = subsample_list[key]
    log(f"Subsample size: {len(sample_ids_sub)}")

    # Subset omics matrices to the shared sample ids
    omics_files = {
        "rna": "rna.tsv",
        "methylation": "methylation.tsv",
        "cnv": "cnv.tsv",
    }

    tmpdir = tempfile.mkdtemp(prefix=f"xintnmf_boot_{iter_id}_", dir=outdir)
    subset_paths = {}

    try:
        for name, fname in omics_files.items():
            full_path = os.path.join(data_dir, fname)
            df = pd.read_csv(full_path, sep="\t", index_col=0)

            missing = sorted(set(sample_ids_sub) - set(df.columns))
            if missing:
                sys.exit(
                    f"Iteration {iter_id}: {len(missing)} sample IDs from the "
                    f"shared subsample are missing from {fname} "
                    f"(e.g. {missing[:3]})"
                )

            # Preserve the shared subsample's sample order across all omics
            df_sub = df[sample_ids_sub]
            sub_path = os.path.join(tmpdir, fname)
            df_sub.to_csv(sub_path, sep="\t")
            subset_paths[name] = sub_path

        # Interaction matrices are feature-feature (e.g. gene x probe,
        # gene x segment), not sample-indexed -- reuse as-is, no subsampling.
        interaction_paths = [
            os.path.join(data_dir, "interaction_rna_methylation.tsv"),
            os.path.join(data_dir, "interaction_rna_cnv.tsv"),
        ]

        os.makedirs(final_out, exist_ok=True)

        cmd = [
            "python", "/opt/X-intNMF/X-intNMF-run.py",
            "--omics_input", subset_paths["rna"], subset_paths["methylation"], subset_paths["cnv"],
            "--interaction_input", *interaction_paths,
            "--output_dir", final_out,
            "--output_format", "csv",
            "--num_components", str(k),
            "--graph_regularization", "1",
            "--max_iter", "5000",
            "--backend", "numpy",
            "--gpu", "-1",
        ]
        log("Running: " + " ".join(cmd))
        subprocess.run(cmd, check=True)

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    log("Done.")


if __name__ == "__main__":
    main()
