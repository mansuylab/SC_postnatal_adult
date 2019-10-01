# Libraries required ----
library(data.table)
library(plgINS)


# Annotation database ----
load("./input/gencode.vM18.anno.RData")


# Importing Salmon files ----
salmon <- salmonImporter(folder = "./input/salmon", anno = anno)


# Calculation of Normalization Factors ----
salmon@norm.factors <- calcNormFactors(salmon@gene.counts)


# PhenoData ----
pheno <- read.delim("input/dataset_pnd_adult.tsv", header = T, stringsAsFactors = F)
rownames(pheno) <- pheno$Samples
pheno <- pheno[colnames(salmon@gene.counts), ]


# Set phenoData of Salmon object ----
salmon@phenoData <- pheno


# Save salmon data as R object ----
save(salmon, file = "./output/SC_controls_rnaseq_salmon.tds.RData", compress = T, compression_level = 3)



# SessionInfo ----
cat("# SessionInfo\n")
devtools::session_info()