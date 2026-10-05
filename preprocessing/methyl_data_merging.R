library(data.table)
library(readxl)

series_matrix_path <- "/mnt/petasan_ccb/alessandro/SCANB/methylation_data/GSE148748_series_matrix.txt"

# Read all lines to locate metadata (sample titles) and the data table 

all_lines <- readLines(series_matrix_path)
# Find the line with sample titles 
title_line_idx <- grep("^!Sample_title", all_lines)
title_line <- all_lines[title_line_idx]

# Split on tabs, remove the field name and surrounding quotes
sample_titles <- strsplit(title_line, "\t")[[1]][-1]
sample_titles <- gsub('^"|"$', "", sample_titles)
sample_titles <- gsub("^Breast cancer ", "", sample_titles)  # -> "PD31029a" etc.

cat("Number of sample titles found:", length(sample_titles), "\n")

# Find the GSM ID order as given in the data table header
table_start_idx <- grep("^!series_matrix_table_begin", all_lines)
header_line <- all_lines[table_start_idx + 1]
gsm_ids <- strsplit(header_line, "\t")[[1]]
gsm_ids <- gsub('^"|"$', "", gsm_ids)
gsm_ids <- gsm_ids[gsm_ids != "ID_REF"]

cat("Number of GSM IDs in table header:", length(gsm_ids), "\n")

# Sanity check: titles and GSM ids should be in the same sample order
stopifnot(length(gsm_ids) == length(sample_titles))
gsm_to_pd <- setNames(sample_titles, gsm_ids)


# Read the beta value table

table_end_idx <- grep("^!series_matrix_table_end", all_lines)
beta_82_raw <- fread(series_matrix_path,
                      sep = "\t",
                      header = TRUE,
                      skip = table_start_idx,        # start right after the marker
                      nrows = (table_end_idx - table_start_idx - 2))

setDF(beta_82_raw)
rownames(beta_82_raw) <- beta_82_raw[[1]]
beta_82_raw[[1]] <- NULL

# Rename columns from GSM IDs to PD-style sample names (column names currently look like "GSM4478173"; map them via gsm_to_pd)
new_names <- gsm_to_pd[colnames(beta_82_raw)]
stopifnot(!any(is.na(new_names)))  # make sure every column was successfully mapped
colnames(beta_82_raw) <- new_names

cat("82-sample matrix dimensions after parsing:", dim(beta_82_raw), "\n")

beta_82_official <- beta_82_raw


# Load the already-normalized 154-sample beta matrix (GSE148906)

path_154_beta <- "/mnt/petasan_ccb/juanra/SCANB/methyl_TNBC/NatComm2025/GSE148906_GPL21145_200408_TNBC154_matrix_normalizedAverageBeta.txt"

dt_154 <- fread(path_154_beta, sep = "\t", header = TRUE)
setDF(dt_154)
rownames(dt_154) <- dt_154[[1]]
dt_154[[1]] <- NULL

beta_154_official <- dt_154

cat("154-sample matrix dimensions:", dim(beta_154_official), "\n")




# Merge datasets

common_cpgs <- intersect(rownames(beta_82_official), rownames(beta_154_official))
cat("Number of CpGs common to both series:", length(common_cpgs), "\n")

beta_82_common  <- beta_82_official[common_cpgs, ]
beta_154_common <- beta_154_official[common_cpgs, ]

# Sanity check: no overlapping sample IDs between the two series
overlap_samples <- intersect(colnames(beta_82_common), colnames(beta_154_common))
if (length(overlap_samples) > 0) {
  cat("WARNING: overlapping sample IDs found between series:", overlap_samples, "\n")
}

merged_beta <- cbind(beta_82_common, beta_154_common)
cat("Merged matrix dimensions (CpGs x samples):", dim(merged_beta), "\n")

# Match with patients in metadata file

#load metadata
ids_match_path <- file.path("/mnt/petasan_ccb/juanra/SCANB/RNAseq/metadata/ids_cohorts_match.xlsx")
ids_match_sheet <- read_excel(ids_match_path, sheet = "1a SCAN-B discovery")
patients <- as.data.frame(ids_match_sheet)
rownames(patients) <- patients[[1]]
patients[[1]] <- NULL
patients <- patients[, 1:4]

merged_beta <- merged_beta[, colnames(merged_beta) %in% rownames(patients)]
cat("Restricted to matched patient set:", dim(merged_beta), "\n")



# Remove CpGs not measured in all samples 
merged_beta <- na.omit(merged_beta)
cat("Final matrix after removing incomplete CpGs:", dim(merged_beta), "\n")

# Save
 
out_path <- "/mnt/petasan_ccb/alessandro/SCANB/methylation_data/matched_methyl_normalized.tsv"
write.table(merged_beta, out_path, sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)
 
