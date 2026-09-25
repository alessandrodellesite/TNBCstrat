library(readxl)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

dt_metadata <- read_excel("/mnt/petasan_ccb/juanra/SCANB/RNAseq/metadata/ids_cohorts_match.xlsx",
                           sheet = "1a SCAN-B discovery")

# clustering solution
rna_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/rna_nmf_clusters.csv")
met_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/methyl_nmf_clusters.csv")
met_results_4 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/methyl_nmf_clusters_4.csv")
cnv_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv/cnv_nmf_clusters.csv")

mofa2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_2.csv")
mofa3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_3.csv")
mofa4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_4.csv")

icluster2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_2.csv")
icluster3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_3.csv")
icluster4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_4.csv")

snf2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_2.csv")
snf3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_3.csv")
snf4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_4.csv")

XintNMF2_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k2_reg.csv")
XintNMF3_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k3_reg.csv")
XintNMF4_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k4_reg.csv")


# Relabel the RAW NMF cluster IDs to the desired display order 
relabel_map <- c("3" = "1", "1" = "2", "2" = "3")
rna_results$Cluster <- factor(relabel_map[as.character(rna_results$Cluster)],
                                levels = c("1", "2", "3"))

#  Organize cluster solutions 
single_omic_list <- list(
  RNA    = rna_results,
  Methylation_k3 = met_results,
  Methylation_k4 = met_results_4,
  CNV    = cnv_results
)

multi_omic_by_k <- list(
  k2 = list(MOFA = mofa2_results, iCluster = icluster2_results, SNF = snf2_results, XintNMF = XintNMF2_results_reg),
  k3 = list(MOFA = mofa3_results, iCluster = icluster3_results, SNF = snf3_results, XintNMF = XintNMF3_results_reg),
  k4 = list(MOFA = mofa4_results, iCluster = icluster4_results, SNF = snf4_results, XintNMF = XintNMF4_results_reg)
)

#  Prepare metadata 
meta <- dt_metadata %>%
  select(PD_ID, OS, OSbin, RFI, RFIbin, DRFI, DRFIbin, Age, TumSize, Grade, LNbinary) %>%
  mutate(across(c(OS, RFI, DRFI, Age, TumSize), as.numeric),
         Grade = as.factor(Grade),
         LNbinary = as.factor(LNbinary))

outcomes <- list(
  OS   = c("OS", "OSbin"),
  RFI  = c("RFI", "RFIbin"),
  DRFI = c("DRFI", "DRFIbin")
)

out_dir  <- "/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/survival"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# function to compute KM + adjusted Cox

run_survival <- function(cluster_df, meta, method_name,
                          time_var, event_var, id_col_cluster = "SampleID") {

  df <- meta %>%
    inner_join(cluster_df, by = c("PD_ID" = id_col_cluster)) %>%
    mutate(Cluster = as.factor(Cluster))

  df_crude <- df %>% filter(!is.na(.data[[time_var]]), !is.na(.data[[event_var]]))
  if (nlevels(droplevels(df_crude$Cluster)) < 2) return(NULL)

  form_crude <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ Cluster"))
  fit <- surv_fit(form_crude, data = df_crude)          # survminer-safe version of survfit()
  logrank_p <- surv_pvalue(fit)$pval

  # adjusted model (drops rows with missing covariates)
  df_adj <- df_crude %>% filter(!is.na(Age), !is.na(TumSize), !is.na(Grade))
  cox_adj_p <- NA_real_
  n_adj <- NA_integer_
  cox_adj <- NULL
  if (nlevels(droplevels(df_adj$Cluster)) >= 2) {
    form_adj <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ Cluster + Age + TumSize + Grade + LNbinary"))
    cox_adj <- tryCatch(coxph(form_adj, data = df_adj), error = function(e) NULL)
    if (!is.null(cox_adj)) {
      cox_adj_p <- anova(cox_adj)["Cluster", "Pr(>|Chi|)"]
      n_adj <- nrow(df_adj)
    }
  }

  list(
    method     = method_name,
    outcome    = time_var,
    n_crude    = nrow(df_crude),
    n_adj      = n_adj,
    n_clusters = nlevels(droplevels(df_crude$Cluster)),
    logrank_p  = logrank_p,
    cox_adj_p  = cox_adj_p,
    fit        = fit,          # kept for later plotting
    df_crude   = df_crude,     # kept for later plotting
    cox_adj    = cox_adj       # kept for later HR extraction if needed
  )
}


# loops over methods x outcomes for one group

run_batch <- function(cluster_list, meta, k_label = NULL) {
  results <- list()
  for (method_name in names(cluster_list)) {
    label <- if (!is.null(k_label)) paste0(method_name, "_", k_label) else method_name
    for (outcome_name in names(outcomes)) {
      time_var  <- outcomes[[outcome_name]][1]
      event_var <- outcomes[[outcome_name]][2]
      key <- paste(label, outcome_name, sep = "_")

      res <- tryCatch(
        run_survival(cluster_list[[method_name]], meta, label, time_var, event_var),
        error = function(e) { message("Failed for ", key, ": ", e$message); NULL }
      )
      results[[key]] <- res
    }
  }
  results
}

# ---- Run for single-omic methods ----
results_single_omic <- run_batch(single_omic_list, meta)

# ---- Run for multi-omic methods, separately per k (fixes the filename collision) ----
results_by_k <- lapply(names(multi_omic_by_k), function(k) {
  run_batch(multi_omic_by_k[[k]], meta, k_label = k)
})
names(results_by_k) <- names(multi_omic_by_k)


# Build summary tables (with BH correction within each group/outcome)

build_summary <- function(results_list) {
  tab <- do.call(rbind, lapply(names(results_list), function(key) {
    res <- results_list[[key]]
    if (is.null(res)) return(NULL)
    data.frame(
      key         = key,
      method      = res$method,
      outcome     = res$outcome,
      n_crude     = res$n_crude,
      n_adj       = res$n_adj,
      n_clusters  = res$n_clusters,
      logrank_p   = res$logrank_p,
      cox_adj_p   = res$cox_adj_p
    )
  }))

  if (is.null(tab) || nrow(tab) == 0) {
    warning("No successful results to summarize.")
    return(tab)
  }

  tab %>%
    group_by(outcome) %>%
    mutate(logrank_p_adj = p.adjust(logrank_p, method = "BH"),
           cox_adj_p_adj = p.adjust(cox_adj_p, method = "BH")) %>%
    ungroup() %>%
    arrange(outcome, cox_adj_p)
}

summary_single_omic <- build_summary(results_single_omic)
summary_k2 <- build_summary(results_by_k$k2)
summary_k3 <- build_summary(results_by_k$k3)
summary_k4 <- build_summary(results_by_k$k4)

summary_all <- bind_rows(
  summary_single_omic %>% mutate(group = "single_omic"),
  summary_k2 %>% mutate(group = "k2"),
  summary_k3 %>% mutate(group = "k3"),
  summary_k4 %>% mutate(group = "k4")
)

# ---- Save summary tables ----
write.csv(summary_single_omic, file.path(out_dir, "summary_single_omic.csv"), row.names = FALSE)
write.csv(summary_k2,          file.path(out_dir, "summary_k2.csv"),          row.names = FALSE)
write.csv(summary_k3,          file.path(out_dir, "summary_k3.csv"),          row.names = FALSE)
write.csv(summary_k4,          file.path(out_dir, "summary_k4.csv"),          row.names = FALSE)
write.csv(summary_all,         file.path(out_dir, "summary_all.csv"),         row.names = FALSE)

# ---- Save full result objects (fits, cox models, data) for later use ----
saveRDS(results_single_omic, file.path(out_dir, "results_single_omic.rds"))
saveRDS(results_by_k,        file.path(out_dir, "results_by_k.rds"))


# Plot only the significant results (based on BH-adjusted p-value)
                        
save_significant_plots <- function(results_list, summary_tab, plot_dir,
                                    p_col = "cox_adj_p_adj", threshold = 0.05) {
  if (is.null(summary_tab) || nrow(summary_tab) == 0) return(invisible(NULL))

  sig_keys <- summary_tab %>%
    filter(!is.na(.data[[p_col]]), .data[[p_col]] < threshold) %>%
    pull(key)

  for (key in sig_keys) {
    res <- results_list[[key]]
    if (is.null(res)) next

    p <- ggsurvplot(res$fit, data = res$df_crude,
                     pval = TRUE, risk.table = TRUE, conf.int = TRUE,
                     xlab = "Years", legend.title = "Cluster",
                     title = paste0(res$method, " — ", res$outcome,
                                    " (adj p=", signif(res$cox_adj_p, 3), ")"))

    fname <- file.path(plot_dir, paste0(key, ".pdf"))
    pdf(fname, width = 7, height = 7)
    print(p)
    dev.off()
  }

  message(length(sig_keys), " significant plot(s) saved from this group.")
}

save_significant_plots(results_single_omic, summary_single_omic, out_dir)
save_significant_plots(results_by_k$k2,      summary_k2,         out_dir)
save_significant_plots(results_by_k$k3,      summary_k3,         out_dir)
save_significant_plots(results_by_k$k4,      summary_k4,         out_dir)




# LAR-referenced survival analysis (TNBCtype4/6 + all k3/k4)

#  Extract Lehman subtypes from metadata 
tnbc4 <- dt_metadata %>%
  select(PD_ID, Cluster = TNBCtype4_n235_notPreCentered) %>%
  filter(!is.na(Cluster), Cluster != "NA")   # catch both true NA and literal string "NA"

tnbc6 <- dt_metadata %>%
  select(PD_ID, Cluster = TNBCtype6_n235_notPreCentered) %>%
  filter(!is.na(Cluster), Cluster != "NA")

cat("Samples in TNBCtype4 clusters:\n")
tnbc4 %>%
  count(Cluster, name = "n_samples") %>%
  print()

cat("\nSamples in TNBCtype6 clusters:\n")
tnbc6 %>%
  count(Cluster, name = "n_samples") %>%
  print()     

#  Group definitions: k3 group and k4 group, each incl. single-omics + TNBCtype 
group_k3 <- list(
  RNA         = list(df = rna_results,          id_col = "SampleID"),
  Methyl_k3      = list(df = met_results,          id_col = "SampleID"),
  CNV         = list(df = cnv_results,          id_col = "SampleID"),
  MOFA_k3     = list(df = mofa3_results,        id_col = "SampleID"),
  iCluster_k3 = list(df = icluster3_results,    id_col = "SampleID"),
  SNF_k3      = list(df = snf3_results,         id_col = "SampleID"),
  XintNMF_k3  = list(df = XintNMF3_results_reg, id_col = "SampleID"),
  TNBCtype4   = list(df = tnbc4,                id_col = "PD_ID"),
  TNBCtype6   = list(df = tnbc6,                id_col = "PD_ID")
)

group_k4 <- list(
  RNA         = list(df = rna_results,          id_col = "SampleID"),
  Methyl_k4      = list(df = met_results_4,          id_col = "SampleID"),
  CNV         = list(df = cnv_results,          id_col = "SampleID"),
  MOFA_k4     = list(df = mofa4_results,        id_col = "SampleID"),
  iCluster_k4 = list(df = icluster4_results,    id_col = "SampleID"),
  SNF_k4      = list(df = snf4_results,         id_col = "SampleID"),
  XintNMF_k4  = list(df = XintNMF4_results_reg, id_col = "SampleID"),
  TNBCtype4   = list(df = tnbc4,                id_col = "PD_ID"),
  TNBCtype6   = list(df = tnbc6,                id_col = "PD_ID")
)

#  Which Cluster label is LAR- matching, per method 
lar_ref_map <- list(
  RNA = "3", Methyl_k3 = "3", Methyl_k4 = "3", CNV = "3",
  MOFA_k3 = "2", MOFA_k4 = "2",
  iCluster_k3 = "2", iCluster_k4 = "4",
  SNF_k3 = "3", SNF_k4 = "2",
  XintNMF_k3 = "2", XintNMF_k4 = "4",
  TNBCtype4 = "LAR", TNBCtype6 = "LAR"
)

#  Function: KM fit + Cox model with LAR forced as reference 
run_survival_LAR <- function(cluster_df, meta, method_name, lar_label,
                              time_var, event_var, id_col_cluster = "SampleID") {

  df <- meta %>%
    inner_join(cluster_df, by = c("PD_ID" = id_col_cluster)) %>%
    mutate(Cluster = as.character(Cluster))

  df_crude <- df %>% filter(!is.na(.data[[time_var]]), !is.na(.data[[event_var]]))
  df_crude$Cluster <- as.factor(df_crude$Cluster)
  if (nlevels(droplevels(df_crude$Cluster)) < 2) return(NULL)
  if (!(lar_label %in% levels(df_crude$Cluster))) {
    message("LAR label '", lar_label, "' not found for ", method_name)
    return(NULL)
  }

  form_crude <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ Cluster"))
  fit <- surv_fit(form_crude, data = df_crude)
  logrank_p <- surv_pvalue(fit)$pval

  df_adj <- df_crude %>% filter(!is.na(Age), !is.na(TumSize), !is.na(Grade))
  df_adj$Cluster <- droplevels(df_adj$Cluster)

  coef_table <- NULL
  n_adj <- NA_integer_
  anova_p <- NA_real_
  if (lar_label %in% levels(df_adj$Cluster) && nlevels(df_adj$Cluster) >= 2) {
    df_adj$Cluster <- relevel(df_adj$Cluster, ref = lar_label)
    form_adj <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ Cluster + Age + TumSize + Grade"))
    cox_adj <- tryCatch(coxph(form_adj, data = df_adj), error = function(e) NULL)
    n_adj <- nrow(df_adj)

    if (!is.null(cox_adj)) {
      s <- summary(cox_adj)
      rows <- grep("^Cluster", rownames(s$coefficients))
      if (length(rows) > 0) {
        coef_table <- data.frame(
          method         = method_name,
          outcome        = time_var,
          lar_ref        = lar_label,
          cluster_vs_lar = sub("^Cluster", "", rownames(s$coefficients)[rows]),
          HR             = s$coefficients[rows, "exp(coef)"],
          lower95        = s$conf.int[rows, "lower .95"],
          upper95        = s$conf.int[rows, "upper .95"],
          p              = s$coefficients[rows, "Pr(>|z|)"]
        )
      }
      # omnibus test for the whole Cluster term (same regardless of reference level)
      anova_p <- tryCatch(anova(cox_adj)["Cluster", "Pr(>|Chi|)"], error = function(e) NA_real_)
    }
  }

  list(
    method     = method_name,
    outcome    = time_var,
    n_crude    = nrow(df_crude),
    n_adj      = n_adj,
    logrank_p  = logrank_p,
    anova_p    = anova_p,
    fit        = fit,
    df_crude   = df_crude,
    coef_table = coef_table
  )
}
                          
# Batch runner for one group (k3 or k4) 
run_group_LAR <- function(group_list, meta, k_name) {
  results <- list()
  for (method_name in names(group_list)) {
    entry <- group_list[[method_name]]
    lar_label <- lar_ref_map[[method_name]]
    for (outcome_name in names(outcomes)) {
      time_var  <- outcomes[[outcome_name]][1]
      event_var <- outcomes[[outcome_name]][2]
      key <- paste(k_name, method_name, outcome_name, sep = "_")

      res <- tryCatch(
        run_survival_LAR(entry$df, meta, method_name, lar_label,
                          time_var, event_var, id_col_cluster = entry$id_col),
        error = function(e) { message("Failed for ", key, ": ", e$message); NULL }
      )
      results[[key]] <- res
    }
  }
  results
}

results_LAR_k3 <- run_group_LAR(group_k3, meta, "k3")
results_LAR_k4 <- run_group_LAR(group_k4, meta, "k4")

#  Build combined p-value table (all methods vs LAR, both groups) 
build_LAR_table <- function(results_list, k_name) {
  tabs <- lapply(names(results_list), function(key) {
    res <- results_list[[key]]
    if (is.null(res) || is.null(res$coef_table)) return(NULL)
    res$coef_table %>%
      mutate(k_group = k_name, n_crude = res$n_crude, n_adj = res$n_adj,
             logrank_p = res$logrank_p, anova_p = res$anova_p)
  })
  bind_rows(tabs)
}

lar_table_k3 <- build_LAR_table(results_LAR_k3, "k3")
lar_table_k4 <- build_LAR_table(results_LAR_k4, "k4")

lar_table_all <- bind_rows(lar_table_k3, lar_table_k4) %>%
  group_by(k_group, outcome) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(k_group, outcome, p)

write.csv(lar_table_all, file.path(out_dir, "LAR_reference_cox_pvalues.csv"), row.names = FALSE)

# Filter to significant comparisons only 
lar_table_significant <- lar_table_all %>%
  filter(p_adj < 0.05) %>%
  arrange(k_group, outcome, p_adj)

lar_table_significant

write.csv(lar_table_significant, file.path(out_dir, "LAR_reference_cox_pvalues_significant.csv"), row.names = FALSE)                        

                        
# Combined KM plots: all methods together, per k and per outcome 
plot_combined_KM <- function(results_list, k_name, outcome_name, out_dir) {
  time_var <- outcomes[[outcome_name]][1]
  keys <- names(results_list)[grepl(paste0("_", time_var, "$"), names(results_list))]

  splots <- list()
  for (key in keys) {
    res <- results_list[[key]]
    if (is.null(res)) next
    splots[[res$method]] <- ggsurvplot(res$fit, data = res$df_crude,
                                        pval = TRUE, conf.int = FALSE,
                                        legend.title = "Cluster",
                                        title = res$method)
  }
  if (length(splots) == 0) return(invisible(NULL))

  combined <- arrange_ggsurvplots(splots, print = FALSE,
                                   ncol = 3, nrow = ceiling(length(splots) / 3))
  fname <- file.path(out_dir, paste0("KM_combined_", k_name, "_", time_var, ".pdf"))
  ggsave(fname, combined, width = 14, height = 4 * ceiling(length(splots) / 3))
}

for (outcome_name in names(outcomes)) {
  plot_combined_KM(results_LAR_k3, "k3", outcome_name, out_dir)
  plot_combined_KM(results_LAR_k4, "k4", outcome_name, out_dir)
}



# forest plots
library(dplyr)
library(ggplot2)

reorder_factor_desc <- function(f, x, fun = median) {
  f <- as.factor(f)
  ord_val <- tapply(x, f, fun, na.rm = TRUE)
  levs <- names(sort(ord_val, decreasing = TRUE))
  factor(f, levels = levs)
}
library(dplyr)
library(ggplot2)

plot_forest_by_k <- function(lar_table, k_name, outcome_name, out_dir) {
  df_plot <- lar_table %>%
    filter(k_group == k_name, outcome == outcome_name)
  
  if (nrow(df_plot) == 0) return(invisible(NULL))
  
  method_order <- c("TNBCtype4", "TNBCtype6", "RNA", "Methyl", "CNV",
                     "MOFA", "SNF", "iCluster", "XintNMF")
  
  df_plot <- df_plot %>%
    mutate(base_method = case_when(
      grepl("^TNBCtype4", method) ~ "TNBCtype4",
      grepl("^TNBCtype6", method) ~ "TNBCtype6",
      grepl("^RNA", method)       ~ "RNA",
      grepl("^Methyl", method)    ~ "Methyl",
      grepl("^CNV", method)       ~ "CNV",
      grepl("^MOFA", method)      ~ "MOFA",
      grepl("^SNF", method)       ~ "SNF",
      grepl("^iCluster", method)  ~ "iCluster",
      grepl("^XintNMF", method)   ~ "XintNMF",
      TRUE ~ method
    )) %>%
    mutate(base_method = factor(base_method, levels = method_order)) %>%
    arrange(base_method, cluster_vs_lar) %>%
    mutate(p_label = paste0("p=", signif(p_adj, 2)),
           facet_label = paste0(method, "  (ANOVA p=", signif(anova_p, 2), ")"))
  
  # lock facet order to match base_method order (not alphabetical)
  facet_order <- df_plot %>% distinct(base_method, facet_label) %>%
    arrange(base_method) %>% pull(facet_label)
  df_plot <- df_plot %>% mutate(facet_label = factor(facet_label, levels = facet_order))
  
  # lock within-facet comparison order (e.g. "1 vs LAR" before "2 vs LAR" before "4 vs LAR")
  # base R equivalent of fct_reorder: build the label, then set factor levels
  # explicitly ordered by the numeric cluster id (descending, so "1" ends up at top of plot)
  df_plot <- df_plot %>%
    mutate(comp_label = paste0(cluster_vs_lar))
  
  ordering <- df_plot %>%
    distinct(facet_label, comp_label, cluster_vs_lar) %>%
    arrange(facet_label, desc(as.numeric(as.character(cluster_vs_lar))))
  
  df_plot <- df_plot %>%
    mutate(comp_label = factor(comp_label, levels = unique(ordering$comp_label)))
  
  p <- ggplot(df_plot, aes(x = HR, y = comp_label)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = p_label, x = upper95), hjust = -0.15, size = 3) +
    facet_grid(rows = vars(facet_label), scales = "free_y", space = "free_y", switch = "y") +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.35))) +
    labs(x = "Hazard Ratio (log scale, vs LAR)", y = NULL,
         title = paste0("Forest plot", k_name, outcome_name),
         caption = "Per-comparison p-values BH-adjusted; ANOVA p = omnibus test for Cluster term") +
    theme_minimal(base_size = 11) +
    theme(
  strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold", margin = margin(0,0,0,0)),
  strip.placement = "outside",
  strip.background = element_blank(),
  panel.spacing = unit(0.15, "lines"),
  axis.text.y = element_text(margin = margin(0,0,0,0)),
  plot.margin = margin(4, 4, 4, 4))
  
  fname <- file.path(out_dir, paste0("forest_", k_name, "_", outcome_name, ".pdf"))
  ggsave(fname, p, width = 9, height = 0.22 * nrow(df_plot) + 0.8 * length(unique(df_plot$facet_label)))
  p
}

for (k_name in c("k3", "k4")) {
  for (outcome_name in names(outcomes)) {
    plot_forest_by_k(lar_table_all, k_name, outcome_name, out_dir)
  }
}


                          
                        
