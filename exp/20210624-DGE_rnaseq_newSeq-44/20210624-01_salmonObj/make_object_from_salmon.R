# Libraries required ----
library(data.table)
library(plgINS)
source("salmonImporter.R")

# Annotation database ----
load("./input/gencode.vM18.anno.RData")


# Importing Salmon files ----
salmon <- salmonImporter1(folder = "./input/salmon", anno = anno)


# Calculation of Normalization Factors ----
salmon@norm.factors <- calcNormFactors(salmon@gene.counts)


# PhenoData ----
p <- read.delim("input/dataset_pnd_adult.tsv", header = T, stringsAsFactors = F)
rownames(p) <- p$SamplesID

pheno <- p
pheno <- pheno[colnames(salmon@gene.counts), ]


# Set phenoData of Salmon object ----
salmon@phenoData <- pheno

# Save salmon data as R object ----
save(salmon, file = "./output/SC_controls_rnaseq_lab_june2021.tds.RData", compress = T, compression_level = 3)


# SessionInfo ----
cat("# SessionInfo\n")
devtools::session_info()
