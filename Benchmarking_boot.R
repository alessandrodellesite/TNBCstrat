#!/usr/bin/env Rscript
## ============================================================
## Full benchmarking driver: runs compute_bootstrap_stability()
## across all method x K combinations and builds a combined,
## K-faceted summary table.
##
## Assumes:
##   - compute_bootstrap_stability() is already defined (see
##     original ALL METHODS TOGETHER chunk)
##   - xintnmf_convert_to_rds.R has already been run, producing
##     xintnmf_bootstrap_k{2,3,4}_clusters/ folders
## ============================================================


library(mclust)
library(pheatmap)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)


# function to compute all metrics

compute_bootstrap_stability <- function(boot_dir, 
                                         reference_clusters,   # named vector: sample_id -> cluster label, full cohort
                                         method_name,
                                         file_pattern = "\\.rds$",
                                         plot = TRUE) {
  
  files <- list.files(boot_dir, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found in ", boot_dir)
  boot_results <- purrr::map(files, readRDS)

  ## integrity checks 
  iter_ids <- purrr::map(boot_results, "iter_id")
  if (any(duplicated(unlist(iter_ids)))) {
    warning(method_name, ": duplicate iter_ids found -- check for resubmitted/duplicate jobs")
  }
  sub_sizes <- purrr::map_int(boot_results, ~ length(.x$sample_ids))
  if (length(unique(sub_sizes)) > 1) {
    warning(method_name, ": inconsistent subsample sizes across replicates -- ",
            paste(names(table(sub_sizes)), table(sub_sizes), sep = "=", collapse = ", "))
  }

  master_samples <- names(reference_clusters) #the full sample cohort from reference_clusters
  n_master <- length(master_samples)

  
  ## Build consensus matrix (considering all bootstrapping replicates)
  
  consensus_sum <- matrix(0, n_master, n_master, dimnames = list(master_samples, master_samples)) #will accumulate how many times 2 samples are clustered together
  connect_count <- matrix(0, n_master, n_master, dimnames = list(master_samples, master_samples)) #how many times did these two samples both appear in the same bootstrap replicate (denominator)
  #loop over every replicate 
  for (res in boot_results) {
    ids  <- res$sample_ids     #replicate id
    labs <- res$cluster       #replicate's cluster assignment
    names(labs) <- ids
    # guard against replicate samples not present in reference_clusters
    ids <- intersect(ids, master_samples)
    labs <- labs[ids]
    same_cluster <- outer(labs, labs, "==") * 1 #builds a matrix comparing every pair of samples' cluster labels (1:same cluster)
    #add it to the initial matrix (only considering rows/cols corresponding to ids)
    consensus_sum[ids, ids] <- consensus_sum[ids, ids] + same_cluster
    connect_count[ids, ids] <- connect_count[ids, ids] + 1
  }
  
  #compute final consensus matrix: for each sample pair, is the fraction of the times they appeared together/ times they land in the same cluster
  consensus_matrix <- consensus_sum / pmax(connect_count, 1)
  diag(consensus_matrix) <- 1

  
  
  ## Cophenetic correlation & dispersion 
  consensus_dist <- as.dist(1 - consensus_matrix) #distance is 1-consensus: high consensus -> low distance
  hc_consensus <- hclust(consensus_dist, method = "average") #hierarchical clustering on that distance
  coph_dist <- cophenetic(hc_consensus) #compute the height of the dendrogram for every pair
  cophenetic_corr  <- cor(consensus_dist, coph_dist, method = "pearson") #Pearson correlation between the actual distances and the tree-implied distances (High: clean clusters)
  dispersion_score <- sum(4 * (consensus_matrix - 0.5)^2) / (n_master^2) #Brunet's dispersion coefficient: value near 1 = consensus values are mostly polarized toward 0 or 1 (confident); near 0 = lots of ambiguous 0.5 values

  ## ARI vs. reference, per replicate 
  ## for each replicate it compares with the reference
  ari_per_rep <- purrr::map_dbl(boot_results, function(res) {
    ids <- intersect(res$sample_ids, master_samples)
    if (length(ids) < 2) return(NA_real_)
    ref_sub <- reference_clusters[ids]
    boot_sub <- res$cluster
    names(boot_sub) <- res$sample_ids
    mclust::adjustedRandIndex(ref_sub, boot_sub[ids]) #computes ARI between the reference partition and this replicate's partition (only considering samples present in the reference)
  })

  ## Jaccard bootstrap stability, per reference cluster (Hennig 2007) ---
  ## For each reference cluster and each bootstrap replicate, find the
  ## best-matching bootstrap cluster (by Jaccard index over member overlap,
  ## restricted to samples present in that replicate) and record that best
  ## Jaccard value. Averaging across replicates gives a per-cluster stability
  ## score -- this is a DIFFERENT question from ARI: ARI asks "does the whole
  ## partition agree", Jaccard-per-cluster asks "does THIS specific cluster
  ## reliably reappear intact". A method can have OK average ARI while one
  ## particular cluster consistently dissolves or merges into others; that
  ## shows up here but is invisible in the ARI/dispersion numbers alone.
  ref_cluster_labels <- sort(unique(reference_clusters))

  #initial matrix
  jaccard_matrix <- matrix(
    NA_real_, nrow = length(boot_results), ncol = length(ref_cluster_labels),
    dimnames = list(NULL, paste0("cluster_", ref_cluster_labels))
  )

  for (rep_i in seq_along(boot_results)) {
    res <- boot_results[[rep_i]]
    ids <- intersect(res$sample_ids, master_samples)
    if (length(ids) < 2) next

    boot_labs <- res$cluster
    names(boot_labs) <- res$sample_ids
    boot_labs <- boot_labs[ids]
    ref_labs  <- reference_clusters[ids]
    boot_cluster_labels <- unique(boot_labs)

    for (cl in ref_cluster_labels) {
      ref_members <- names(ref_labs)[ref_labs == cl]
      if (length(ref_members) == 0) next  # this reference cluster absent from this replicate's ids

      best_jaccard <- 0
      for (bcl in boot_cluster_labels) {
        boot_members <- names(boot_labs)[boot_labs == bcl]
        inter <- length(intersect(ref_members, boot_members))
        union_n <- length(union(ref_members, boot_members))
        jacc <- if (union_n == 0) 0 else inter / union_n
        if (jacc > best_jaccard) best_jaccard <- jacc
      }
      jaccard_matrix[rep_i, paste0("cluster_", cl)] <- best_jaccard
    }
  }

  jaccard_per_cluster <- colMeans(jaccard_matrix, na.rm = TRUE)
  jaccard_sd_per_cluster <- apply(jaccard_matrix, 2, sd, na.rm = TRUE)
  mean_jaccard <- mean(jaccard_per_cluster, na.rm = TRUE)  # overall, unweighted across clusters

  ## --- Per-sample stability ---
  per_sample_stability <- purrr::map_dbl(master_samples, function(s) {
    own_cluster <- reference_clusters[s]
    same_members <- names(reference_clusters)[reference_clusters == own_cluster & names(reference_clusters) != s]
    same_members <- intersect(same_members, master_samples)
    mean(consensus_matrix[s, same_members], na.rm = TRUE)
  })
  names(per_sample_stability) <- master_samples

  ## --- Optional heatmap ---
  if (plot) {
    ord <- order(reference_clusters)
    pheatmap(consensus_matrix[ord, ord],
             cluster_rows = FALSE, cluster_cols = FALSE,
             annotation_row = data.frame(Cluster = factor(reference_clusters[ord]),
                                          row.names = master_samples[ord]),
             show_rownames = FALSE, show_colnames = FALSE,
             color = colorRampPalette(c("white", "#1A73E8"))(50),
             main = sprintf("%s consensus matrix (%d replicates)", method_name, length(boot_results)))
  }

  list(
    method              = method_name,
    n_replicates        = length(boot_results),
    consensus_matrix    = consensus_matrix,
    cophenetic_corr     = cophenetic_corr,
    dispersion_score    = dispersion_score,
    ari_per_rep         = ari_per_rep,
    mean_ari            = mean(ari_per_rep, na.rm = TRUE),
    median_ari          = median(ari_per_rep, na.rm = TRUE),
    sd_ari              = sd(ari_per_rep, na.rm = TRUE),
    jaccard_matrix      = jaccard_matrix,
    jaccard_per_cluster = jaccard_per_cluster,
    jaccard_sd_per_cluster = jaccard_sd_per_cluster,
    mean_jaccard        = mean_jaccard,
    per_sample_stability = per_sample_stability
  )
}

boot_base <- "/mnt/petasan_ccb/alessandro/SCANB/multiomics"

## --- Helper: load a reference CSV into a named vector (sample_id -> cluster) ---
load_reference <- function(path) {
  df <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)
  stopifnot(all(c("SampleID", "Cluster") %in% names(df)))
  ref <- df$Cluster
  names(ref) <- df$SampleID
  ref
}

## --- Config grid: one row per method x K combination that actually exists ---
config <- tribble(
  ~method,    ~k, ~boot_dir,                                                                  ~ref_file,                                                              ~file_pattern,
  "iCluster+", 2,  file.path(boot_base, "bootstrap", "icluster", "icluster_bootstrap_k2"),      file.path(boot_base, "output_icluster", "iclusters_clusters_2.csv"),   "\\.rds$",
  "iCluster+", 3,  file.path(boot_base, "bootstrap", "icluster", "icluster_bootstrap_k3"),      file.path(boot_base, "output_icluster", "iclusters_clusters_3.csv"),   "\\.rds$",
  "iCluster+", 4,  file.path(boot_base, "bootstrap", "icluster", "icluster_bootstrap_k4"),      file.path(boot_base, "output_icluster", "iclusters_clusters_4.csv"),   "\\.rds$",
  "MOFA",     2,  file.path(boot_base, "bootstrap", "mofa", "mofa_bootstrap_k2"),              file.path(boot_base, "output_mofa", "mofa_km_clusters_2.csv"),        "\\.rds$",
  "MOFA",     3,  file.path(boot_base, "bootstrap", "mofa", "mofa_bootstrap_k3"),              file.path(boot_base, "output_mofa", "mofa_km_clusters_3.csv"),        "\\.rds$",
  "MOFA",     4,  file.path(boot_base, "bootstrap", "mofa", "mofa_bootstrap_k4"),              file.path(boot_base, "output_mofa", "mofa_km_clusters_4.csv"),        "\\.rds$",
  "SNF",      2,  file.path(boot_base, "bootstrap", "snf", "snf_bootstrap_k2"),                file.path(boot_base, "output_snf", "snf_clusters_2.csv"),             "\\.rds$",
  "SNF",      3,  file.path(boot_base, "bootstrap", "snf", "snf_bootstrap_k3"),                file.path(boot_base, "output_snf", "snf_clusters_3.csv"),             "\\.rds$",
  "SNF",      4,  file.path(boot_base, "bootstrap", "snf", "snf_bootstrap_k4"),                file.path(boot_base, "output_snf", "snf_clusters_4.csv"),             "\\.rds$",
  "X-intNMF",  2,  file.path(boot_base, "bootstrap", "xintnmf", "xintnmf_bootstrap_k2_clusters"), file.path(boot_base, "output_xintnmf", "xintNMF_clusters_k2_reg.csv"), "\\.rds$",
  "X-intNMF",  3,  file.path(boot_base, "bootstrap", "xintnmf", "xintnmf_bootstrap_k3_clusters"), file.path(boot_base, "output_xintnmf", "xintNMF_clusters_k3_reg.csv"), "\\.rds$",
  "X-intNMF",  4,  file.path(boot_base, "bootstrap", "xintnmf", "xintnmf_bootstrap_k4_clusters"), file.path(boot_base, "output_xintnmf", "xintNMF_clusters_k4_reg.csv"), "\\.rds$"
)


## --- Run stability analysis for every row in the config grid ---
consensus_dir <- "/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/bootstrapping/consensus_matrices"
dir.create(consensus_dir, showWarnings = FALSE, recursive = TRUE)

all_stability <- purrr::pmap(config, function(method, k, boot_dir, ref_file, file_pattern) {
  reference_clusters <- load_reference(ref_file)
  message(sprintf("Running %s k=%d ...", method, k))

  png_path <- file.path(consensus_dir, sprintf("consensus_%s_k%d.png", method, k))
  png(png_path, width = 1200, height = 1000, res = 150)

  res <- compute_bootstrap_stability(
    boot_dir           = boot_dir,
    reference_clusters = reference_clusters,
    method_name        = sprintf("%s (k=%d)", method, k),
    file_pattern       = file_pattern,
    plot               = TRUE
  )

  dev.off()   

  res$K <- k
  res$method_base <- method
  res
})





## --- Cross-method replicate-alignment sanity check ---
## Confirms that, for a given K, all methods' bootstrap files were built
## from the SAME shared subsamples (as intended) -- i.e. that the
## *sample sets* per iter_id match across methods, since that's what
## the shared boot_subsamples.rds is supposed to guarantee.

check_shared_subsamples <- function(config_subset) {
  file_lists <- purrr::pmap(config_subset, function(method, k, boot_dir, ref_file, file_pattern) {
    files <- list.files(boot_dir, pattern = file_pattern, full.names = TRUE)
    res <- purrr::map(files, readRDS)
    purrr::set_names(purrr::map(res, ~ sort(.x$sample_ids)), purrr::map_chr(res, ~ as.character(.x$iter_id)))
  })
  names(file_lists) <- config_subset$method

  common_iters <- purrr::reduce(purrr::map(file_lists, names), intersect)
  mismatches <- 0
  for (it in common_iters) {
    sets <- purrr::map(file_lists, ~ .x[[it]])
    sets <- purrr::compact(sets)
    if (length(unique(sets)) > 1) mismatches <- mismatches + 1
  }
  message(sprintf("  Checked %d shared iter_ids: %d mismatched sample sets across methods",
                   length(common_iters), mismatches))
  invisible(mismatches)
}

for (k_val in c(2, 3, 4)) {
  message(sprintf("Cross-method subsample check for K=%d:", k_val))
  check_shared_subsamples(config[config$k == k_val, ])
}

## --- Combined summary table, faceted by K ---
summary_table <- purrr::map_dfr(all_stability, function(s) {
  data.frame(
    Method            = s$method_base,
    K                 = s$K,
    N_Replicates      = s$n_replicates,
    Mean_ARI          = round(s$mean_ari, 3),
    Median_ARI        = round(s$median_ari, 3),
    SD_ARI            = round(s$sd_ari, 3),
    #Mean_Jaccard      = round(s$mean_jaccard, 3),
    Cophenetic        = round(s$cophenetic_corr, 3),
    Dispersion        = round(s$dispersion_score, 3)
  )
})

summary_table <- summary_table[order(summary_table$K, summary_table$Method), ]
print(summary_table)

## Optional: write out for reporting
write.csv(summary_table, file.path(boot_base, "bootstrap_stability_summary.csv"),  row.names = FALSE)



## PLOTS


# Boxplots of the full per-replicate ARI distribution 

ari_long <- purrr::map_dfr(all_stability, function(s) {
  data.frame(Method = s$method_base, K = s$K, ARI = s$ari_per_rep)
})
p_box_ari <- ggplot(ari_long, aes(x = Method, y = ARI, fill = Method)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~ K, labeller = labeller(K = function(x) paste0("K = ", x))) +
  labs(title = "Distribution of per-replicate ARI vs reference",
       y = "ARI", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/bootstrapping",
               "boxplots_ari.png"),
    width = 1400, height = 1000, res = 150)
print(p_box_ari)
dev.off()                                         
                                      

# Boxplots of per-cluster Jaccard values (pooled across clusters) 
## Each point here is one (bootstrap replicate x reference cluster) best-match
## Jaccard value, so this shows the spread Hennig's method is meant to reveal:
## whether instability is concentrated in one or two clusters or spread evenly.
                                      
jaccard_long <- purrr::map_dfr(all_stability, function(s) {
  jm <- s$jaccard_matrix
  df <- as.data.frame(jm)
  df$Method <- s$method_base
  df$K <- s$K
  tidyr::pivot_longer(df, cols = starts_with("cluster_"),
                       names_to = "RefCluster", values_to = "Jaccard")
})
p_box_jaccard <- ggplot(jaccard_long, aes(x = Method, y = Jaccard, fill = Method)) +
  geom_boxplot(outlier.size = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = c(0.6, 0.85), linetype = "dashed", color = "grey40") +
  facet_wrap(~ K, labeller = labeller(K = function(x) paste0("K = ", x))) +
  labs(title = "Distribution of per-cluster Jaccard stability",
       y = "Best-match Jaccard", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
print(p_box_jaccard)

png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/bootstrapping",
               "boxplots_jaccard.png"),
    width = 1400, height = 1000, res = 150)
print(p_box_jaccard)
dev.off()                                        

## Per-cluster Jaccard breakdown (which SPECIFIC clusters are unstable) ---
## This is the plot that actually answers "which cluster is the problem".
## IMPORTANT: cluster labels (cluster_1, cluster_2, ...) are arbitrary,
## method-specific indices with NO cross-method correspondence -- iCluster's
## "cluster_1" and SNF's "cluster_1" are not the same biological group. So
## methods are put in SEPARATE facet panels (facet_grid) rather than dodged
## side-by-side at the same x-position, which would visually (and wrongly)
## imply the labels line up across methods.
                                      
# not useful to make comparisons across methods, more to see within each method which clusters are most conserved                                     
                                      
p_cluster_detail <- ggplot(jaccard_long, aes(x = RefCluster, y = Jaccard, fill = RefCluster)) +
  geom_boxplot(outlier.size = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = c(0.6, 0.85), linetype = "dashed", color = "grey40") +
  facet_grid(K ~ Method, scales = "free_x", space = "free_x",
             labeller = labeller(K = function(x) paste0("K = ", x))) +
  labs(title = "Per-cluster Jaccard stability, by method and reference cluster",
       y = "Best-match Jaccard") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cluster_detail)
       

png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/bootstrapping",
               "boxplots_jaccard_percluster.png"),
    width = 1400, height = 1000, res = 150)
print(p_cluster_detail)
dev.off()    


# Visualizing cophenetic correlation, dispersion, and per-sample stability                                 

## Cophenetic correlation & dispersion score, by method x K 
## These summarize consensus-matrix "cleanliness" (how tree-like/polarized
## the pairwise co-clustering structure is) but were previously only in
## the summary table as numbers. Plotting them alongside ARI/Jaccard makes
## it possible to spot cases where these three lines of evidence disagree
## (e.g. high ARI but low dispersion -- consensus is polarized but the
## dendrogram doesn't reflect it cleanly).
                                 
diag_long <- summary_table %>%
  select(Method, K, Cophenetic_Corr, Dispersion) %>%
  tidyr::pivot_longer(cols = c(Cophenetic_Corr, Dispersion),
                       names_to = "Metric", values_to = "Value")


p_diag <- ggplot(diag_long, aes(x = Method, y = Value, fill = Method)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = round(Value, 2)), vjust = -0.4, size = 3.2) +
  facet_grid(Metric ~ K, labeller = labeller(K = function(x) paste0("K = ", x))) +
  ylim(0, 1.08) +
  labs(title = "Cophenetic correlation & dispersion scores by method and k",
       y = NULL, x = NULL) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
print(p_diag)

png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/bootstrapping",
               "cophenetic_dispersion.png"),
    width = 1400, height = 1000, res = 150)
print(p_diag)
dev.off()
                                             
                                             
                                             
# QUESTO NON SO SE HA SENSO
                                             
## --- 7. Per-sample stability: full distribution, by method x K ---
per_sample_long <- purrr::map_dfr(all_stability, function(s) {
  data.frame(Method = s$method_base, K = s$K,
             SampleID = names(s$per_sample_stability),
             Stability = s$per_sample_stability)
})

p_per_sample <- ggplot(per_sample_long, aes(x = Method, y = Stability, fill = Method)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  facet_wrap(~ K, labeller = labeller(K = function(x) paste0("K = ", x))) +
  ylim(0, 1) +
  labs(title = "Per-sample stability distribution",
       subtitle = "Dashed line at 0.5: samples below this are ambiguous/boundary cases",
       y = "Mean consensus with own cluster", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
print(p_per_sample)

png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/bootstrapping",
               "persample_stability.png"),
    width = 1400, height = 1000, res = 150)
print(p_per_sample)
dev.off()  
                                      

## --- 8. Summary: fraction of "unstable" samples per method x K ---
## Complements the violin plot with a single interpretable number:
## how many samples are borderline/unreliable cluster members.
unstable_summary <- per_sample_long %>%
  group_by(Method, K) %>%
  summarise(
    n_samples          = n(),
    frac_below_0.5     = mean(Stability < 0.5, na.rm = TRUE),
    frac_below_0.7     = mean(Stability < 0.7, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(K, Method)
print(unstable_summary)


# QUESTO SOTTO NON MI PIACE MOLTO 


#4 — Paired test between methods (exploiting the shared-subsample design)
# Since ari_per_rep inside compute_bootstrap_stability isn't tagged with iter_id, pairing purely by vector 
# position risks silently misaligning replicates if file listing order differs between methods. I wrote a small 
# standalone extractor that re-derives ARI values explicitly keyed by iter_id, so the pairing is guaranteed correct 
# regardless of file ordering:

## --- Extract ARI per replicate, named by iter_id (for safe pairing) ---
get_named_ari <- function(boot_dir, reference_clusters, file_pattern = "\\.rds$") {
  files <- list.files(boot_dir, pattern = file_pattern, full.names = TRUE)
  boot_results <- purrr::map(files, readRDS)
  master_samples <- names(reference_clusters)

  ari_vec <- purrr::map_dbl(boot_results, function(res) {
    ids <- intersect(res$sample_ids, master_samples)
    if (length(ids) < 2) return(NA_real_)
    ref_sub  <- reference_clusters[ids]
    boot_sub <- res$cluster
    names(boot_sub) <- res$sample_ids
    mclust::adjustedRandIndex(ref_sub, boot_sub[ids])
  })
  names(ari_vec) <- purrr::map_chr(boot_results, ~ as.character(.x$iter_id))
  ari_vec
}

## --- Build named-ARI vectors for every method x K combination ---
ari_by_config <- purrr::pmap(config, function(method, k, boot_dir, ref_file, file_pattern) {
  reference_clusters <- load_reference(ref_file)
  list(method = method, K = k,
       ari = get_named_ari(boot_dir, reference_clusters, file_pattern))
})

## --- Pairwise paired Wilcoxon signed-rank test, within each K ---
## Valid because every method was run on the SAME 1000 subsamples (verified
## earlier by check_shared_subsamples), so ARI values at matching iter_id
## are naturally paired -- a paired test is more powerful here than an
## unpaired comparison of the two ARI distributions.
pairwise_results <- list()

for (k_val in unique(config$k)) {
  methods_here <- purrr::keep(ari_by_config, ~ .x$K == k_val)
  if (length(methods_here) < 2) next

  combos <- combn(seq_along(methods_here), 2, simplify = FALSE)
  for (pair in combos) {
    a <- methods_here[[pair[1]]]
    b <- methods_here[[pair[2]]]

    common_iters <- intersect(names(a$ari), names(b$ari))
    valid <- !is.na(a$ari[common_iters]) & !is.na(b$ari[common_iters])
    common_iters <- common_iters[valid]
    if (length(common_iters) < 10) next   # not enough paired replicates to test

    x <- a$ari[common_iters]
    y <- b$ari[common_iters]
    wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE))

    pairwise_results[[length(pairwise_results) + 1]] <- data.frame(
      K            = k_val,
      Method_A     = a$method,
      Method_B     = b$method,
      n_paired     = length(common_iters),
      median_ARI_A = median(x),
      median_ARI_B = median(y),
      median_diff  = median(x - y),
      W_statistic  = unname(wt$statistic),
      p_value      = wt$p.value
    )
  }
}

pairwise_df <- dplyr::bind_rows(pairwise_results)

## Multiple-testing correction, applied within each K (each K is treated
## as its own family of comparisons rather than pooling corrections globally)
pairwise_df <- pairwise_df %>%
  group_by(K) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(K, p_adj)

print(pairwise_df)

## --- Visualize: heatmap of adjusted p-values per K ---
p_pairwise <- ggplot(pairwise_df, aes(x = Method_A, y = Method_B, fill = p_adj)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif(p_adj, 2)), size = 3) +
  scale_fill_gradient(low = "#D62828", high = "#F0F0F0", limits = c(0, 1),
                       name = "BH-adj. p") +
  facet_wrap(~ K, labeller = labeller(K = function(x) paste0("K = ", x))) +
  labs(title = "Paired Wilcoxon test: pairwise ARI stability comparisons",
       subtitle = "Red = significant difference in stability between methods (BH-adjusted p < 0.05)",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_pairwise)
                                      
png(file.path("/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/bootstrapping",
               "pairwise_test.png"),
    width = 1400, height = 1000, res = 150)
print(p_pairwise)
dev.off()                                      
