library(readxl)
library(dplyr)
library(survival)
library(survminer)

# load file with survival scores
dt_metadata <- read_excel("/mnt/petasan_ccb/juanra/SCANB/RNAseq/metadata/ids_cohorts_match.xlsx", sheet = "1a SCAN-B discovery")

#load clustering solutions
rna_results   <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_rna/rna_nmf_clusters.csv")
met_results   <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_meth/methyl_nmf_clusters.csv")
cnv_results   <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/singleomic/nmf_cnv/cnv_nmf_clusters.csv")

#mofa
mofa2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_2.csv")
mofa3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_3.csv")
mofa4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_mofa/mofa_km_clusters_4.csv")

#icluster
icluster2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_2.csv")
icluster3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_3.csv")
icluster4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_icluster/iclusters_clusters_4.csv")

#snf
snf2_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_2.csv")
snf3_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_3.csv")
snf4_results <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_snf/snf_clusters_4.csv")

#xintnmf
XintNMF2_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k2_reg.csv")
XintNMF3_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k3_reg.csv")
XintNMF4_results_reg <- read.csv("/mnt/petasan_ccb/alessandro/SCANB/multiomics/output_xintnmf/xintNMF_clusters_k4_reg.csv")


single_omic_list <- list(
  RNA    = rna_results,
  Methyl = met_results,
  CNV    = cnv_results
)

multi_omic_by_k <- list(
  k2 = list(
    MOFA    = mofa2_results,
    iCluster = icluster2_results,
    SNF     = snf2_results,
    XintNMF = XintNMF2_results_reg
  ),
  k3 = list(
    MOFA    = mofa3_results,
    iCluster = icluster3_results,
    SNF     = snf3_results,
    XintNMF = XintNMF3_results_reg
  ),
  k4 = list(
    MOFA    = mofa4_results,
    iCluster = icluster4_results,
    SNF     = snf4_results,
    XintNMF = XintNMF4_results_reg
  )
)

meta <- dt_metadata %>%
  select(PD_ID, OS, OSbin, RFI, RFIbin, DRFI, DRFIbin, Age, TumSize, Grade) %>%
  mutate(across(c(OS, RFI, DRFI, Age, TumSize), as.numeric),
         Grade = as.factor(Grade))

# function that computes all survivals

plot_dir <- "/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/survival"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

run_survival <- function(cluster_df, meta, method_name,
                          time_var, event_var, id_col_cluster = "SampleID") {
  
  df <- meta %>%
    inner_join(cluster_df, by = c("PD_ID" = id_col_cluster)) %>%
    mutate(Cluster = as.factor(Cluster))
  
  df_crude <- df %>% filter(!is.na(.data[[time_var]]), !is.na(.data[[event_var]]))
  if (nlevels(droplevels(df_crude$Cluster)) < 2) return(NULL)
  
  form_crude <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ Cluster"))
  fit  <- survfit(form_crude, data = df_crude)
  cox_crude <- coxph(form_crude, data = df_crude)
  logrank_p <- surv_pvalue(fit)$pval
  
  # adjusted model (drops rows with missing covariates)
  df_adj <- df_crude %>% filter(!is.na(Age), !is.na(TumSize), !is.na(Grade))
  cox_adj_p <- NA_real_
  n_adj <- NA_integer_
  if (nlevels(droplevels(df_adj$Cluster)) >= 2) {
    form_adj <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ Cluster + Age + TumSize + Grade"))
    cox_adj <- tryCatch(coxph(form_adj, data = df_adj), error = function(e) NULL)
    if (!is.null(cox_adj)) {
      # overall p-value for the Cluster term (likelihood ratio test, all cluster levels)
      cox_adj_p <- anova(cox_adj)["Cluster", "Pr(>|Chi|)"]
      n_adj <- nrow(df_adj)
    }
  }
  
  p <- ggsurvplot(fit, data = df_crude,
                   pval = TRUE,
                   risk.table = TRUE,
                   conf.int = TRUE,
                   xlab = "Years",
                   legend.title = "Cluster",
                   title = paste0(method_name, " — ", time_var))
  
  # save plot
  fname <- file.path(plot_dir, paste0(method_name, "_", time_var, ".pdf"))
  pdf(fname, width = 7, height = 7)
  print(p)
  dev.off()
  
  list(
    method = method_name,
    outcome = time_var,
    n_crude = nrow(df_crude),
    n_adj = n_adj,
    n_clusters = nlevels(droplevels(df_crude$Cluster)),
    logrank_p = logrank_p,
    cox_adj_p = cox_adj_p,
    cox_crude_summary = summary(cox_crude),
    plot = p
  )
}


outcomes <- list(
  OS   = c("OS", "OSbin"),
  RFI  = c("RFI", "RFIbin"),
  DRFI = c("DRFI", "DRFIbin")
)

run_batch <- function(cluster_list, meta) {
  results <- list()
  for (method_name in names(cluster_list)) {
    for (outcome_name in names(outcomes)) {
      time_var  <- outcomes[[outcome_name]][1]
      event_var <- outcomes[[outcome_name]][2]
      key <- paste(method_name, outcome_name, sep = "_")
      
      res <- tryCatch(
        run_survival(cluster_list[[method_name]], meta, method_name, time_var, event_var),
        error = function(e) {
          message("Failed for ", key, ": ", e$message)
          NULL
        }
      )
      results[[key]] <- res
    }
  }
  results
}

results_single_omic <- run_batch(single_omic_list, meta)

results_by_k <- lapply(multi_omic_by_k, run_batch, meta = meta)
# results_by_k$k2, results_by_k$k3, results_by_k$k4



build_summary <- function(results_list) {
  tab <- do.call(rbind, lapply(names(results_list), function(key) {
    res <- results_list[[key]]
    if (is.null(res)) return(NULL)
    data.frame(
      method = res$method,
      outcome = res$outcome,
      n_crude = res$n_crude,
      n_adj = res$n_adj,
      n_clusters = res$n_clusters,
      logrank_p = res$logrank_p,
      cox_adj_p = res$cox_adj_p
    )
  }))
  
  tab %>%
    group_by(outcome) %>%
    mutate(logrank_p_adj = p.adjust(logrank_p, method = "BH"),
           cox_adj_p_adj = p.adjust(cox_adj_p, method = "BH")) %>%
    ungroup() %>%
    arrange(outcome, logrank_p)
}

summary_single_omic <- build_summary(results_single_omic)

summary_k2 <- build_summary(results_by_k$k2)
summary_k3 <- build_summary(results_by_k$k3)
summary_k4 <- build_summary(results_by_k$k4)           
                        
out_dir <- "/mnt/petasan_ccb/alessandro/SCANB/plots/benchmarking/survival"

write.csv(summary_single_omic, file.path(out_dir, "summary_single_omic.csv"), row.names = FALSE)
write.csv(summary_k2,          file.path(out_dir, "summary_k2.csv"),          row.names = FALSE)
write.csv(summary_k3,          file.path(out_dir, "summary_k3.csv"),          row.names = FALSE)
write.csv(summary_k4,          file.path(out_dir, "summary_k4.csv"),          row.names = FALSE)


summary_all <- bind_rows(
  summary_single_omic %>% mutate(group = "single_omic"),
  summary_k2 %>% mutate(group = "k2"),
  summary_k3 %>% mutate(group = "k3"),
  summary_k4 %>% mutate(group = "k4")
)

write.csv(summary_all, file.path(out_dir, "summary_all.csv"), row.names = FALSE)

saveRDS(results_single_omic, file.path(out_dir, "results_single_omic.rds"))
saveRDS(results_by_k,        file.path(out_dir, "results_by_k.rds"))                        
