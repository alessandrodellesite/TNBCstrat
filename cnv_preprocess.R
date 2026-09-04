library(biomaRt)
library(GenomicRanges)

# Data loading and matching samples with rnaseq dataset

# copy number variants file
cnv_data <- read.csv("/mnt/petasan_ccb/juanra/SCANB/matched_data/matched_cnv.tsv", sep = "\t", header = TRUE) 
summary(cnv_data$Seg.CN)

#metadata for ploidy values
library(readxl)
metadata <- read_excel("/mnt/petasan_ccb/juanra/SCANB/RNAseq/metadata/ids_cohorts_match.xlsx", sheet = "1a SCAN-B discovery")
metadata <- as.data.frame(metadata)
ploidy_lookup <- setNames(as.numeric(metadata$ASCAT_PLOIDY), metadata$PD_ID) #ploidy values x sample


# rnaseq file 
rna_data <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.rds")

dim(cnv_data) #64084 x 6
dim(rna_data) #15268 x 235

common_samples <- intersect(colnames(rna_data), unique(cnv_data$Sample))
length(common_samples)
rna_data <- rna_data[, common_samples]

# Retrieve coordinates (start/end pos, chr) of each gene in rnaseq matrix - USE GRCH37
#mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://grch37.ensembl.org")
#mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)

#Retry loop
mart <- NULL
attempts <- 0
while (is.null(mart) && attempts < 5) {
  attempts <- attempts + 1
  mart <- tryCatch({
    useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)
  }, error = function(e) {
    message("Attempt ", attempts, " failed: ", conditionMessage(e))
    Sys.sleep(5)
    NULL
  })
}

#get attributes: Ensembl ID, chr, start, end
gene_coords <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", 
                                    "start_position", "end_position"),
                     filters = "ensembl_gene_id",
                     values = rownames(rna_data), 
                     mart = mart)
dim(gene_coords) # it was able to find the coordinates of 15115 genes againts 15268

# Converted gene_coords into GenomicRanges object
gene_gr <- GRanges(seqnames = gene_coords$chromosome_name,
                   ranges = IRanges(start = gene_coords$start_position, 
                                    end = gene_coords$end_position),
                   gene_id = gene_coords$ensembl_gene_id)

#EMPTY matrix for DNA (genes x samples)
dna_matrix <- matrix(NA, nrow = nrow(gene_coords), ncol = length(common_samples),
                     dimnames = list(gene_coords$ensembl_gene_id, common_samples))


# associate segments to genes that are within its genomic coordinates:

# matrix to keep track of how many segments fall within one gene's coordinates
n_segments_matrix <- matrix(0L,
                             nrow = nrow(gene_coords),
                             ncol = length(common_samples),
                             dimnames = list(gene_coords$ensembl_gene_id, common_samples))

# matrices to preserve min/max copy number per gene per sample (GDC-style extra info)
min_cn_matrix <- matrix(NA_real_, nrow = nrow(gene_coords), ncol = length(common_samples),
                         dimnames = list(gene_coords$ensembl_gene_id, common_samples))
max_cn_matrix <- matrix(NA_real_, nrow = nrow(gene_coords), ncol = length(common_samples),
                         dimnames = list(gene_coords$ensembl_gene_id, common_samples))

# collector for discordant (gene, sample) records - built up across the whole loop
discordant_records <- list()

#for loop that iterates over all patient and for each segment it re-computes the total allele count and finds the genes present in its coordinates 
for (samp in common_samples) {
  # Subset CNV for this sample
  samp_cnv <- cnv_data[cnv_data$Sample == samp, ]
  
    # PLOIDY CORRECTION
  
  #samp_ploidy <- ploidy_lookup[[samp]]
  #stopifnot(!is.na(samp_ploidy))   # fail if a sample is missing ploidy
  ## ploidy-adjust Seg.CN BEFORE back-transforming 
  #samp_cnv$Seg.CN_adj <- samp_cnv$Seg.CN - log2(samp_ploidy / 2)
  #samp_cnv$alleles <- round(2 * (2^samp_cnv$Seg.CN_adj))   #this has to replace the line below

  # log2 Seg.CN converted back to alleles (total copy number)
  samp_cnv$alleles <- round(2 * (2^samp_cnv$Seg.CN))

  # Map segments to genes using Genomic Overlap
  seg_gr <- GRanges(seqnames = samp_cnv$Chromosome,
                    ranges = IRanges(start = samp_cnv$Start.Position,
                                     end = samp_cnv$End.Position),
                    alleles = samp_cnv$alleles)

  # Find which segment/s are present within the gene's coordinates
  hits <- findOverlaps(gene_gr, seg_gr)

  # computes intersection between the segment-gene pairs (how much every segment covers the gene)
  intersections <- pintersect(gene_gr[queryHits(hits)], seg_gr[subjectHits(hits)])
  overlap_widths <- width(intersections)

  # temporary dataframe: one row per (gene, segment) overlap pair
  hit_df <- data.frame(
    gene_idx = queryHits(hits),
    alleles  = mcols(seg_gr)$alleles[subjectHits(hits)],
    overlap  = overlap_widths
  )

  # QC: count how many segments overlap each gene for this sample
  seg_counts <- aggregate(overlap ~ gene_idx, data = hit_df, FUN = length)
  n_segments_matrix[seg_counts$gene_idx, samp] <- seg_counts$overlap

  # identify multi-segment genes, and among those, which are genuinely discordant
  multi_idx <- seg_counts$gene_idx[seg_counts$overlap > 1]
  discordant_idx <- integer(0)

  if (length(multi_idx) > 0) {
    hit_df_multi <- hit_df[hit_df$gene_idx %in% multi_idx, ]
    allele_variance <- aggregate(alleles ~ gene_idx, data = hit_df_multi, FUN = var)
    discordant_idx <- allele_variance$gene_idx[allele_variance$alleles > 0]

    if (length(discordant_idx) > 0) {
      # keep a full record of these discordant events for later review
      rec <- hit_df_multi[hit_df_multi$gene_idx %in% discordant_idx, ]
      rec$gene_id <- gene_coords$ensembl_gene_id[rec$gene_idx]
      rec$sample  <- samp
      gene_total_width <- width(gene_gr)[rec$gene_idx]
      rec$pct_of_gene  <- round(100 * rec$overlap / gene_total_width, 1)
      discordant_records[[length(discordant_records) + 1]] <- rec
    }
  }

  # weighted average, applied to EVERY gene (1 segment, concordant, or discordant)
  # weighted mean of allele values, weighted by how much of the gene each segment covers
  weighted_vals <- hit_df |>
    dplyr::group_by(gene_idx) |>
    dplyr::summarise(
      weighted_cn = round(sum(alleles * overlap) / sum(overlap)),
      min_cn      = min(alleles),
      max_cn      = max(alleles),
      .groups = "drop"
    )

  dna_matrix[weighted_vals$gene_idx, samp]   <- weighted_vals$weighted_cn
  min_cn_matrix[weighted_vals$gene_idx, samp] <- weighted_vals$min_cn
  max_cn_matrix[weighted_vals$gene_idx, samp] <- weighted_vals$max_cn
}

dim(dna_matrix)
sum(is.na(dna_matrix))   # should now only be genes with ZERO overlapping segments

#  Discordance summary tables: keeps track of the genes that had more that one segment within its coordinates

discordant_all <- do.call(rbind, discordant_records)

if (!is.null(discordant_all)) {

  # per (gene, sample): what did discordance look like for this specific patient
  per_sample_discordance <- discordant_all |>
    dplyr::group_by(gene_id, sample) |>
    dplyr::summarise(
      n_segments_this_sample = dplyr::n(),
      values_this_sample     = paste(sort(unique(alleles)), collapse = ","),
      n_distinct_this_sample = dplyr::n_distinct(alleles),
      min_pct_this_sample    = min(pct_of_gene),
      max_pct_this_sample    = max(pct_of_gene),
      .groups = "drop"
    )

  # gene-level: aggregated FROM the per-sample table, not from raw pooled rows
  discordant_summary <- per_sample_discordance |>
    dplyr::group_by(gene_id) |>
    dplyr::summarise(
      n_samples_discordant = dplyr::n(),
      n_distinct_patterns  = dplyr::n_distinct(values_this_sample),
      example_patterns     = paste(unique(values_this_sample)[1:min(3, dplyr::n_distinct(values_this_sample))], collapse = " | "),
      mean_min_pct         = round(mean(min_pct_this_sample), 1),
      mean_max_pct         = round(mean(max_pct_this_sample), 1),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_samples_discordant))

} else {
  per_sample_discordance <- data.frame()
  discordant_summary <- data.frame()
}

head(discordant_summary, 20)
nrow(discordant_summary)


#Per-gene: in how many samples did it overlap >1 segment?
n_multi_per_gene <- rowSums(n_segments_matrix > 1)
# per ogni gene: # samples che hanno almeno 1 caso di più di un segmento per quel gene

#summary table
multi_seg_summary <- data.frame(
  gene_id  = names(n_multi_per_gene),
  n_samples = n_multi_per_gene
) |> dplyr::arrange(desc(n_samples))

cat("Genes overlapping >1 segment in at least 1 sample:",
    sum(n_multi_per_gene > 0), "\n")
cat("Genes overlapping >1 segment in >10% of samples:",
    sum(n_multi_per_gene > 0.1 * length(common_samples)), "\n")

head(multi_seg_summary, 20)

# How many NAs were introduced by discordance
cat("NA values in dna_matrix:", sum(is.na(dna_matrix)), "\n")
cat("Proportion NA:", round(mean(is.na(dna_matrix)), 4), "\n")


#extreme_genes <- rownames(dna_matrix)[apply(dna_matrix, 1, function(x) any(x > 200, na.rm = TRUE))]
#extreme_genes
#VOPP1, EGFR, and LANCL2 are known to be amplified in breast cancer, particularly as a coordinated cluster on the short arm of chromosome 7 (7p11.2) - therefore not excluded

# BINARY DNA MATRICES

amp_mat    <- dna_matrix > 5 # T means that gene in that sample is amplifier (more than 5 alleles)
homdel_mat <- dna_matrix == 0 # T means that gene in that sample is deleted (0 copies)
hetdel_mat <- dna_matrix == 1 # T means that gene in that sample has 1 copy

# gets total true values
n_true_vals <- function(mat){
  true_vals <- which(mat == TRUE, arr.ind = TRUE)
  return(length(true_vals))
}

n_true_vals(amp_mat)
n_true_vals(homdel_mat)
n_true_vals(hetdel_mat)

# BINARY RNA MATRICES
#fit distribution for each gene and find outliers
#T: gene in that samples is in 5 or 95% quantile of expression 


# All false matrices (genes x samples)
overexp_mat  <- matrix(FALSE, nrow = nrow(rna_data), ncol = ncol(rna_data), 
                       dimnames = list(rownames(rna_data), colnames(rna_data)))
underexp_mat <- matrix(FALSE, nrow = nrow(rna_data), ncol = ncol(rna_data), 
                       dimnames = list(rownames(rna_data), colnames(rna_data)))

# genes found in both DNA and RNA matrices
common_genes <- intersect(rownames(dna_matrix), rownames(rna_data)) #15226

library(MASS)
for (gene in common_genes) {
  exp_values <- as.numeric(rna_data[gene, ])
  # Fit distribution and get 5%/95% thresholds
  est <- fitdistr(na.omit(exp_values), "normal")$estimate
  qvalues <- qnorm(c(0.05, 0.95), est[1], est[2])
    
  overexp_mat[gene, ]  <- exp_values > qvalues[2]
  underexp_mat[gene, ] <- exp_values < qvalues[1]
} 
#n_true_vals(overexp_mat)
#n_true_vals(underexp_mat)

# Match dimensions
overexp_mat  <- overexp_mat[common_genes, common_samples]
underexp_mat <- underexp_mat[common_genes, common_samples]
amp_mat      <- amp_mat[common_genes, common_samples]
homdel_mat   <- homdel_mat[common_genes, common_samples]
hetdel_mat   <- hetdel_mat[common_genes, common_samples]


#final outlier Matrices DA RIVEDERE VALORI
amp_over_95quantile <- amp_mat & overexp_mat #132526
del_under_5quantile <- homdel_mat & underexp_mat #666
het_under_5quantile <- hetdel_mat & underexp_mat #72662

library(tibble)
final_counts <- tibble(
  Gene = common_genes,
  Amp_Over_95quantile = rowSums(amp_over_95quantile, na.rm = TRUE),
  Del_Under_5quantile = rowSums(del_under_5quantile, na.rm = TRUE),
  het_under_5quantile = rowSums(het_under_5quantile, na.rm= TRUE)
)


# Choose minimum number of samples for a gene that have to have an amplification + overexpression for feature selection

#Dataframe with TRUE genes retained applying different thresholds for each binary dataset

thresholds <- 1:10

#counts for each row
conteggi_amp <- rowSums(amp_over_95quantile, na.rm = TRUE)
conteggi_del <- rowSums(del_under_5quantile, na.rm = TRUE)
conteggi_het <- rowSums(het_under_5quantile, na.rm = TRUE)

#vectors
riga_amp <- sapply(thresholds, function(t) sum(conteggi_amp >= t, na.rm = TRUE))
riga_del <- sapply(thresholds, function(t) sum(conteggi_del >= t, na.rm = TRUE))
riga_het <- sapply(thresholds, function(t) sum(conteggi_het >= t, na.rm = TRUE))

valori <- rbind(riga_amp, riga_del, riga_het)
rownames(valori) <- NULL
valori_con_totale <- rbind(valori, colSums(valori))

df_finale <- data.frame(
  Threshold = c("amp_over_95quantile", "del_under_5quantile", "het_under_5quantile", "Total"),
  valori_con_totale
)
colnames(df_finale)[2:11] <- thresholds

print(df_finale, row.names = FALSE)


#Threshold: >7 samples for amp_over_95quantile and het_under_5quantile

# genes from each group with same threshold for samples
thr <- 7
cont_amp <- rowSums(amp_over_95quantile, na.rm = TRUE)
geni_amp <- names(cont_amp[cont_amp > thr])

cont_del <- rowSums(del_under_5quantile, na.rm = TRUE)
geni_del <- names(cont_del[cont_del > 0]) # homozygous deletions are rarer events I don't want to under-select

cont_het <- rowSums(het_under_5quantile, na.rm = TRUE)
geni_het <- names(cont_het[cont_het > thr])

# gather unique genes
geni_final <- union(union(geni_amp, geni_del), geni_het)

# check
cat("Genes from AMP:", length(geni_amp), "\n")
cat("Genes from DEL:", length(geni_del), "\n")
cat("Genes from HET:", length(geni_het), "\n")
cat("Total genes combined:", length(geni_final), "\n")

# Filter matrices with final genes
dna_matrix_filt <- dna_matrix[geni_final, , drop = FALSE]
rna_data_filt   <- rna_data[geni_final, , drop = FALSE]
rna_matrix_filt <- as.matrix(rna_data_filt)

#genes with at least 1 NA
sum(rowSums(is.na(dna_matrix_filt)) > 0) 

# Calcola la media di ogni riga ESCLUDENDO gli NA (altrimenti la media sarebbe NA)
# Se una riga fosse completamente vuota, la media sarebbe NaN; usiamo un'accortezza.


sum(is.na(dna_matrix_filt))                  # total missing cells
mean(is.na(dna_matrix_filt))                 # proportion missing across the whole matrix
table(rowSums(is.na(dna_matrix_filt)))       # distribution: how many genes have 1 NA, 2 NA, ... all NA?

                   
row_means <- rowMeans(dna_matrix_filt, na.rm = TRUE)

# Trova le coordinate (riga, colonna) di tutti gli NA nella matrice
na_idx <- which(is.na(dna_matrix_filt), arr.ind = TRUE)

# Sostituisci ogni NA con la media della sua rispettiva riga
# na_idx[, 1] estrae l'indice di riga di ciascun valore mancante
dna_matrix_filt[na_idx] <- row_means[na_idx[, 1]]

#save
saveRDS(dna_matrix_filt, "/mnt/petasan_ccb/alessandro/SCANB/cnv_processed.rds")
