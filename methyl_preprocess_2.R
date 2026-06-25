# Methylation data preprocessing

library(flexmix)
library(parallel)
library(lattice)
library(readxl)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ELMER)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(GenomeInfoDb)
library(GenomicRanges)

# Point to the root mount where the data is stored
adjusted_data <- readRDS("/mnt/petasan_ccb/alessandro/SCANB/adjusted_data.rds")
exp.data <- read.table("/mnt/petasan_ccb/alessandro/SCANB/rna_logtransformed.tsv", header=TRUE, sep="\t", row.names=1)
                                       
# Filtering to exclude chrX/Y localization and non-CpG probes

print("Filtering probes")
options(sesameData.offline = TRUE)
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

keep_autosomes <- !(ann$chr %in% c("chrX", "chrY"))
keep_cpg <- grepl("^cg", ann$Name)
probes_to_keep <- ann$Name[keep_autosomes & keep_cpg]
probes_filtered <- adjusted_data[rownames(adjusted_data) %in% probes_to_keep, ]
print("Probes filtered")

# ELMER

print("ELMER analysis started")

met_filtered <- as.matrix(probes_filtered)
exp_matrix <- as.matrix(exp.data)
clean_gene_names <- gsub("\\..*", "", rownames(exp_matrix))
rownames(exp_matrix) <- clean_gene_names

common_samples <- intersect(colnames(met_filtered), colnames(exp_matrix))
met_filtered   <- met_filtered[, common_samples]
exp_matrix     <- exp_matrix[, common_samples]

sample_meta <- data.frame(
  group = rep("Cancer", length(common_samples)),
  row.names = common_samples,
  stringsAsFactors = FALSE
)

print("Get TSS coordinates")
gene_coords <- getTSS(genome = "hg19")
# take the TSS from our total gene list
match_idx <- match(rownames(exp_matrix), gene_coords$ensembl_gene_id)
valid_genes <- !is.na(match_idx)


distal_probe_coords <- get.feature.probe(genome = "hg19", met.platform = "EPIC")
# take the distal probes from our total CpGs list
final_distal_probes <- intersect(rownames(met_filtered), names(distal_probe_coords))


#exploit complete annotation already done
all_probes_gr <- GRanges(
  seqnames = ann$chr,
  ranges = IRanges(start = ann$pos, end = ann$pos),
  strand = ann$strand,
  Name = ann$Name
)
names(all_probes_gr) <- ann$Name

# promoters (+/- 2kb dal TSS)
promoter_regions <- promoters(gene_coords, upstream = 2000, downstream = 2000)

# take the probes from our total CpGs list that fall into this region
overlaps <- findOverlaps(all_probes_gr, promoter_regions)
probes_promoter_names <- unique(names(all_probes_gr[queryHits(overlaps)]))
final_promoter_list <- intersect(rownames(met_filtered), probes_promoter_names)


print("Create the multiassaysexperiments")
# genes summarizedExperiment
se_exp <- SummarizedExperiment(
  assays = list(exp = exp_matrix[valid_genes, ]),
  rowRanges = gene_coords[match_idx[valid_genes],]
)

seqlevelsStyle(se_exp) <- "UCSC"

# ENHANCERS Multiassay experiment
se_met_en <- SummarizedExperiment(
  assays = list(met = met_filtered[final_distal_probes, ]),
  rowRanges = all_probes_gr[final_distal_probes]
)

seqlevelsStyle(se_met_en) <- "UCSC"

mae_enhancer <- MultiAssayExperiment(
  experiments = list("DNA methylation" = se_met_en, "Gene expression" = se_exp),
  colData = sample_meta
)


# PROMOTERS Multiassay experiment
se_met_pr <- SummarizedExperiment(
  assays = list(met = met_filtered[final_promoter_list, ]),
  rowRanges = all_probes_gr[final_promoter_list]
)
seqlevelsStyle(se_met_pr) <- "UCSC"

mae_promoter <- MultiAssayExperiment(
  experiments = list("DNA methylation" = se_met_pr, "Gene expression" = se_exp),
  colData = sample_meta
)


# script to run get.pair() with 1000 permutation and filtering 50k probes for variance

# Parameters
n_probes <- 50000  
n_permu  <- 1000
n_cores  <- 30      


# ENHANCERS 

print("Enhancers ELMER analysis started")
                                       
probe_sd_en <- apply(assay(mae_enhancer, "DNA methylation"), 1, sd)
top_probes_en <- names(sort(probe_sd_en, decreasing = TRUE))[1:n_probes]

subset_list_en <- list("DNA methylation" = top_probes_en, 
                       "Gene expression" = rownames(mae_enhancer[["Gene expression"]]))

mae_enhancer_top <- mae_enhancer[subset_list_en, , ]

nearGenes_enhancer <- GetNearGenes(data = mae_enhancer_top, 
                                   probes = top_probes_en, 
                                   numFlankingGenes = 20)
print("Enhancers get.pair analysis started")
                                       
pairs_enhancer <- get.pair(
  data = mae_enhancer,
  nearGenes = nearGenes_enhancer,
  group.col = "group",
  group1 = "Cancer",
  group2 = "Cancer",
  mode = "unsupervised", 
  diff.dir = "both",
  permu.size = n_permu,      
  filter.probes = FALSE, 
  cores = n_cores
)

saveRDS(pairs_enhancer, "result_pairs_enhancer.rds")

print("Enhancers get.pair analysis finished")

# PROMOTERS 
print("Promoter ELMER analysis started")

probe_sd_pr <- apply(assay(mae_promoter, "DNA methylation"), 1, sd)
top_probes_pr <- names(sort(probe_sd_pr, decreasing = TRUE))[1:n_probes]

# Create new MAE filtering for these
subset_list_pr <- list("DNA methylation" = top_probes_pr, 
                       "Gene expression" = rownames(mae_promoter[["Gene expression"]]))

mae_promoter_top <- mae_promoter[subset_list_pr, , ]

nearGenes_promoters <- GetNearGenes(data = mae_promoter_top, 
                                    probes = top_probes_pr,
                                    numFlankingGenes = 2 
)
                                       
print("Promoters get.pair analysis started")

#get.pair (per le correlazioni)
pairs_promoter <- get.pair(
  data = mae_promoter,
  nearGenes = nearGenes_promoters,
  group.col = "group",
  group1 = "Cancer",
  group2 = "Cancer",
  mode = "unsupervised", 
  diff.dir = "hypo",
  permu.size = n_permu,
  filter.probes = FALSE, 
  cores = n_cores
) 

saveRDS(pairs_promoter, "result_pairs_promoters.rds")

print("Promoters get.pair analysis finished!")
print("All done!")
