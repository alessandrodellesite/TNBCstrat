library(readxl)
library(dplyr)
library(survival)
library(survminer)


# load file with survival scores
dt_metadata <- read_excel("ids_cohorts_match.xlsx", sheet = "1a SCAN-B discovery") 

#load clustering solutions
rna_clusters   <- read.csv("rnaseq_nmf_clusters.csv")
met_clusters   <- read.csv("rnaseq_nmf_clusters.csv")
cnv_clusters   <- read.csv("rnaseq_nmf_clusters.csv")
mofa_clusters  <- read.csv("rnaseq_nmf_clusters.csv")
snf_clusters   <- read.csv("rnaseq_nmf_clusters.csv")
icl_clusters   <- read.csv("rnaseq_nmf_clusters.csv")
xint_clusters  <- read.csv("rnaseq_nmf_clusters.csv")






