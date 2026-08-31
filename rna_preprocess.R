# RNAseq data loading and preprocessing

dt <- read.csv("/mnt/petasan_ccb/juanra/SCANB/matched_data/matched_rna.tsv", sep="\t", row.names = 1)
#sum(dt < 0) #no neg vals

# FPKM to TPM
tpm_matrix <- apply(dt, 2, function(x) x / sum(as.numeric(x))*10^6)

#keep genes with a sum of expression across samples > 0 
data_filtered <- tpm_matrix[rowSums(dt) > 0, ]

# FILTERING FOR LOW EXPRESSION DATA 
threshold_tpm <- 2

# min number of samples 
min_samples <- 0.05*dim(data_filtered)[2]

#keep genes where at least 5% of samples have a TPM > 2
keep <- rowSums(data_filtered > threshold_tpm) >= min_samples

# Apply the filter to TPM data
data_fil <- data_filtered[keep, ]

# log2 transformation
data_prepped <- log2(data_fil + 1)



#Visualization plots before-after filtering
library(ggplot2)
library(patchwork) 

df_before <- data.frame(mean_expr = rowMeans(log2(data_filtered + 1)))
df_after <- data.frame(mean_expr = rowMeans(log2(data_fil + 1)))

plot_before <- ggplot(df_before, aes(x = mean_expr)) +
  geom_density(fill = "firebrick", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Before Filtering", x = "Log2(TPM + 1)", y = "Density") +
  geom_vline(xintercept = log2(threshold_tpm + 1), linetype = "dashed")

plot_after <- ggplot(df_after, aes(x = mean_expr)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  theme_minimal() +
  labs(title = "After Filtering", x = "Log2(TPM + 1)", y = "Density") +
  geom_vline(xintercept = log2(threshold_tpm + 1), linetype = "dashed")

#plot_before + plot_after
#save plots
combined_plot <- plot_before + plot_after
ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/expr_before_after_filtering.png",
       plot = combined_plot,
       width = 10, height = 5, dpi = 300)
                    
                    
# Function to get mean and sd for plotting
get_stats <- function(df, label) {
  log_data <- log2(df + 1)
  data.frame(
    mean = rowMeans(log_data),
    sd = apply(log_data, 1, sd),
    status = label
  )
  }


stats_before <- get_stats(data_filtered, "Before Filtering")
stats_after  <- get_stats(data_fil, "After Filtering")

plot_mv_before <- ggplot(stats_before, aes(x = mean, y = sd)) +
  geom_point(alpha = 0.1, color = "firebrick") +
  geom_smooth(method = "gam", color = "black") + # Adds a trend line
  theme_minimal() +
  labs(title = "Mean-Variance: Before", x = "Mean Log2(TPM+1)", y = "Std. Deviation")

plot_mv_after <- ggplot(stats_after, aes(x = mean, y = sd)) +
  geom_point(alpha = 0.1, color = "steelblue") +
  geom_smooth(method = "gam", color = "black") +
  theme_minimal() +
  labs(title = "Mean-Variance: After", x = "Mean Log2(TPM+1)", y = "Std. Deviation")

#plot_mv_before + plot_mv_after
                    
#save plots
combined_plot_mv <- plot_mv_before + plot_mv_after
ggsave("/mnt/petasan_ccb/alessandro/SCANB/plots/mv_before_after_filtering.png",
       plot = combined_plot_mv,
       width = 10, height = 5, dpi = 300)
                    
#dim(tpm_matrix)
#dim(data_prepped)

#clean ensamble nomenclature
head(rownames(data_prepped))
clean_ids <- gsub("\\..*", "", rownames(data_prepped))
rownames(data_prepped) <- clean_ids

#save
saveRDS(data_prepped, "/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.rds")
