library(readxl)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

dt_metadata <- read_excel("/mnt/petasan_ccb/juanra/SCANB/RNAseq/metadata/ids_cohorts_match.xlsx",
                           sheet = "1a SCAN-B discovery")

rna_results   <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/rna_nmf_clusters.csv")
met_results   <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/methyl_nmf_clusters.csv")
met_results_4 <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/methyl_nmf_clusters_4.csv")
cnv_results   <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv/cnv_nmf_clusters.csv")

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

#  Relabel RNA clusters to desired display order 
relabel_map <- c("3" = "1", "1" = "2", "2" = "3")
rna_results$Cluster <- factor(relabel_map[as.character(rna_results$Cluster)], levels = c("1", "2", "3"))

# Lehman subtypes from metadata 
tnbc4 <- dt_metadata %>%
  select(PD_ID, Cluster = TNBCtype4_n235_notPreCentered) %>%
  filter(!is.na(Cluster), Cluster != "NA")

tnbc6 <- dt_metadata %>%
  select(PD_ID, Cluster = TNBCtype6_n235_notPreCentered) %>%
  filter(!is.na(Cluster), Cluster != "NA")

#metadata 
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

out_dir <- "/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/survival"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# One group definition per k, each entry: df, id_col, lar_label 
groups <- list(
  k2 = list(
    MOFA     = list(df = mofa2_results,        id_col = "SampleID", lar_label = NULL),
    iCluster = list(df = icluster2_results,    id_col = "SampleID", lar_label = NULL),
    SNF      = list(df = snf2_results,         id_col = "SampleID", lar_label = NULL),
    XintNMF  = list(df = XintNMF2_results_reg, id_col = "SampleID", lar_label = NULL)
  ),
  k3 = list(
    RNA         = list(df = rna_results,          id_col = "SampleID", lar_label = "3"),
    Methyl_k3   = list(df = met_results,          id_col = "SampleID", lar_label = "3"),
    CNV         = list(df = cnv_results,          id_col = "SampleID", lar_label = "3"),
    MOFA_k3     = list(df = mofa3_results,        id_col = "SampleID", lar_label = "2"),
    iCluster_k3 = list(df = icluster3_results,    id_col = "SampleID", lar_label = "2"),
    SNF_k3      = list(df = snf3_results,         id_col = "SampleID", lar_label = "3"),
    XintNMF_k3  = list(df = XintNMF3_results_reg, id_col = "SampleID", lar_label = "2"),
    TNBCtype4   = list(df = tnbc4,                id_col = "PD_ID",    lar_label = "LAR"),
    TNBCtype6   = list(df = tnbc6,                id_col = "PD_ID",    lar_label = "LAR")
  ),
  k4 = list(
    RNA         = list(df = rna_results,          id_col = "SampleID", lar_label = "3"),
    Methyl_k4   = list(df = met_results_4,        id_col = "SampleID", lar_label = "3"),
    CNV         = list(df = cnv_results,          id_col = "SampleID", lar_label = "3"),
    MOFA_k4     = list(df = mofa4_results,        id_col = "SampleID", lar_label = "2"),
    iCluster_k4 = list(df = icluster4_results,    id_col = "SampleID", lar_label = "4"),
    SNF_k4      = list(df = snf4_results,         id_col = "SampleID", lar_label = "2"),
    XintNMF_k4  = list(df = XintNMF4_results_reg, id_col = "SampleID", lar_label = "4"),
    TNBCtype4   = list(df = tnbc4,                id_col = "PD_ID",    lar_label = "LAR"),
    TNBCtype6   = list(df = tnbc6,                id_col = "PD_ID",    lar_label = "LAR")
  )
)

# Single function: KM + logrank + omnibus adjusted Cox ANOVA + pairwise Cox (reference: LAR)

run_survival_full <- function(cluster_df, meta, method_name, lar_label,
                               time_var, event_var, id_col_cluster = "SampleID") {

  df <- meta %>%
    inner_join(cluster_df, by = c("PD_ID" = id_col_cluster)) %>%
    mutate(Cluster = as.character(Cluster))

  df_crude <- df %>% filter(!is.na(.data[[time_var]]), !is.na(.data[[event_var]]))
  df_crude$Cluster <- as.factor(df_crude$Cluster)
  if (nlevels(droplevels(df_crude$Cluster)) < 2) return(NULL)

  form_crude <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ Cluster"))
  fit <- surv_fit(form_crude, data = df_crude)
  logrank_p <- surv_pvalue(fit)$pval

  df_adj <- df_crude %>% filter(!is.na(Age), !is.na(TumSize), !is.na(Grade), !is.na(LNbinary))
  df_adj$Cluster <- droplevels(df_adj$Cluster)

  n_adj <- NA_integer_
  anova_p <- NA_real_
  coef_table <- NULL

  if (nlevels(df_adj$Cluster) >= 2) {
    n_adj <- nrow(df_adj)
    form_adj <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ Cluster + Age + TumSize + Grade + LNbinary"))

    #  Omnibus ANOVA (no reference )
    cox_default <- tryCatch(coxph(form_adj, data = df_adj), error = function(e) NULL)
    if (!is.null(cox_default)) {
      anova_p <- tryCatch(anova(cox_default)["Cluster", "Pr(>|Chi|)"], error = function(e) NA_real_)
    }

    # LAR-referenced pairwise coefficients, only if a lar_label was supplied (so for k=3 and 4)
    if (!is.null(lar_label) && lar_label %in% levels(df_adj$Cluster)) {
      df_lar <- df_adj
      df_lar$Cluster <- relevel(df_lar$Cluster, ref = lar_label)
      cox_lar <- tryCatch(coxph(form_adj, data = df_lar), error = function(e) NULL)
      if (!is.null(cox_lar)) {
        s <- summary(cox_lar)
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
      }
    }
  }

  list(
    method     = method_name,
    outcome    = time_var,
    n_crude    = nrow(df_crude),
    n_adj      = n_adj,
    n_clusters = nlevels(droplevels(df_crude$Cluster)),
    logrank_p  = logrank_p,
    anova_p    = anova_p,
    fit        = fit,
    df_crude   = df_crude,
    coef_table = coef_table
  )
}


# Batch runner for one k-group
run_group <- function(group_list, meta, k_name) {
  results <- list()
  for (method_name in names(group_list)) {
    entry <- group_list[[method_name]]
    for (outcome_name in names(outcomes)) {
      time_var  <- outcomes[[outcome_name]][1]
      event_var <- outcomes[[outcome_name]][2]
      key <- paste(k_name, method_name, outcome_name, sep = "_")

      res <- tryCatch(
        run_survival_full(entry$df, meta, method_name, entry$lar_label,
                           time_var, event_var, id_col_cluster = entry$id_col),
        error = function(e) { message("Failed for ", key, ": ", e$message); NULL }
      )
      results[[key]] <- res
    }
  }
  results
}

results_by_k <- lapply(names(groups), function(k) run_group(groups[[k]], meta, k))
names(results_by_k) <- names(groups)


# Table 1: ANOVA omnibus p-value per clustering solution (no correction)
                       
build_anova_table <- function(results_list, k_name) {
  do.call(rbind, lapply(names(results_list), function(key) {
    res <- results_list[[key]]
    if (is.null(res)) return(NULL)
    data.frame(
      k_group    = k_name,
      method     = res$method,
      outcome    = res$outcome,
      n_crude    = res$n_crude,
      n_adj      = res$n_adj,
      n_clusters = res$n_clusters,
      logrank_p  = res$logrank_p,
      anova_p    = res$anova_p
    )
  }))
}

anova_table_all <- bind_rows(lapply(names(results_by_k), function(k) build_anova_table(results_by_k[[k]], k))) %>%
  arrange(k_group, outcome, anova_p)

write.csv(anova_table_all, file.path(out_dir, "anova_table_all.csv"), row.names = FALSE)

anova_significant <- anova_table_all %>% filter(!is.na(anova_p), anova_p < 0.05)

cat("\nSignificant ANOVA results (raw p < 0.05, no correction) \n")
print(anova_significant %>% select(k_group, method, outcome, n_adj, anova_p))

write.csv(anova_significant, file.path(out_dir, "anova_table_significant.csv"), row.names = FALSE)

                                    
# Table 2: LAR vs all other clusters, per clustering solution

build_lar_table <- function(results_list, k_name) {
  tabs <- lapply(names(results_list), function(key) {
    res <- results_list[[key]]
    if (is.null(res) || is.null(res$coef_table)) return(NULL)
    res$coef_table %>%
      mutate(k_group = k_name, n_crude = res$n_crude, n_adj = res$n_adj,
             logrank_p = res$logrank_p, anova_p = res$anova_p)
  })
  bind_rows(tabs)
}

lar_table_all <- bind_rows(lapply(names(results_by_k), function(k) build_lar_table(results_by_k[[k]], k))) %>%
  arrange(k_group, outcome, method, p)

write.csv(lar_table_all, file.path(out_dir, "LAR_reference_cox_pvalues.csv"), row.names = FALSE)


                                  
# Forest plots (raw p, raw anova_p)
                                  
plot_forest_by_k <- function(lar_table, k_name, outcome_name, out_dir) {
  df_plot <- lar_table %>% filter(k_group == k_name, outcome == outcome_name)
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
    mutate(p_label = paste0("p=", signif(p, 2)),
           facet_label = paste0(method, "  (ANOVA p=", signif(anova_p, 2), ")"))

  facet_order <- df_plot %>% distinct(base_method, facet_label) %>%
    arrange(base_method) %>% pull(facet_label)
  df_plot <- df_plot %>% mutate(facet_label = factor(facet_label, levels = facet_order))

  df_plot <- df_plot %>% mutate(comp_label = cluster_vs_lar)
  ordering <- df_plot %>%
    distinct(facet_label, comp_label, cluster_vs_lar) %>%
    arrange(facet_label, desc(as.numeric(as.character(cluster_vs_lar))))
  df_plot <- df_plot %>% mutate(comp_label = factor(comp_label, levels = unique(ordering$comp_label)))

  p <- ggplot(df_plot, aes(x = HR, y = comp_label)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = p_label, x = upper95), hjust = -0.15, size = 3) +
    facet_grid(rows = vars(facet_label), scales = "free_y", space = "free_y", switch = "y") +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.35))) +
    labs(x = "Hazard Ratio (vs LAR)", y = NULL,
         title = paste0("Forest plot — ", k_name, " — ", outcome_name),
         caption = "p: covariate-adjusted Cox model; ANOVA p: omnibus test across clusters") +
    theme_minimal(base_size = 11) +
    theme(
      strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold", margin = margin(0,0,0,0)),
      strip.placement = "outside",
      strip.background = element_blank(),
      panel.spacing = unit(0.15, "lines"),
      axis.text.y = element_text(margin = margin(0,0,0,0)),
      plot.margin = margin(4, 4, 4, 4)
    )

  fname <- file.path(out_dir, paste0("forest_", k_name, "_", outcome_name, ".pdf"))
  ggsave(fname, p, width = 9, height = 0.22 * nrow(df_plot) + 0.8 * length(unique(df_plot$facet_label)))
  p
}

for (k_name in c("k3", "k4")) {
  for (outcome_name in names(outcomes)) {
    plot_forest_by_k(lar_table_all, k_name, outcome_name, out_dir)
  }
}


# Combined KM plots
plot_combined_KM <- function(results_list, k_name, outcome_name, out_dir) {
  time_var <- outcomes[[outcome_name]][1]
  keys <- names(results_list)[grepl(paste0("_", time_var, "$"), names(results_list))]

  splots <- list()
  for (key in keys) {
    res <- results_list[[key]]
    if (is.null(res)) next
    splots[[res$method]] <- {
      n_cl <- length(res$fit$strata)
      labs <- sub(".*=", "", names(res$fit$strata))   # "cluster=1" -> "1"
      
      p <- ggsurvplot(
        res$fit, data = res$df_crude,
        pval = TRUE, conf.int = FALSE,
        legend.title = "",            # removes the cluster title
        legend.labs = paste("Cluster", labs),   # or just `labs` for maximum space saving
        legend = "top",
        title = res$method
      )
      
      p$plot <- p$plot +
      guides(colour = guide_legend(nrow = if (n_cl >= 4) 2 else 1, byrow = TRUE)) +
        theme(
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"),
          legend.spacing.x = unit(0.1, "cm")
        )
      p
                                   }
    
  }
  if (length(splots) == 0) return(invisible(NULL))

  combined <- arrange_ggsurvplots(splots, print = FALSE,
                                   ncol = 3, nrow = ceiling(length(splots) / 3))
  fname <- file.path(out_dir, paste0("KM_combined_", k_name, "_", time_var, ".pdf"))
  ggsave(fname, combined, width = 14, height = 4 * ceiling(length(splots) / 3))
}

for (k_name in names(results_by_k)) {
  for (outcome_name in names(outcomes)) {
    plot_combined_KM(results_by_k[[k_name]], k_name, outcome_name, out_dir)
  }
}

#save
saveRDS(results_by_k, file.path(out_dir, "results_by_k.rds"))
