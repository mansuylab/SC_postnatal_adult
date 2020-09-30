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
p1 <- read.delim("input/metadata.txt", header = T, stringsAsFactors = F)
p1 <- p1[p1$Assay.Type == "RS", c(2,8)]
colnames(p1) <- c("Group", "Samples")
p1$Samples <- gsub(pattern = "_1", replacement = "", x = p1$Samples)
rownames(p1) <- p1$Samples
p1 <- p1[p1$Group %in% c("PND14", "PNW8"),]

p2 <- read.delim("input/dataset_pnd_adult.tsv", header = T, stringsAsFactors = F)
rownames(p2) <- p2$Samples
p2 <- p2[p2$Group != "Adult",]

pheno <- rbind(p2, p1)
pheno <- pheno[colnames(salmon@gene.counts), ]


# Set phenoData of Salmon object ----
salmon@phenoData <- pheno
salmon@phenoData$Group1 <- "PND14_PND15"
salmon@phenoData$Group1[salmon@phenoData$Group == "PNW8"] <- "PNW8"
salmon@phenoData$Group1[salmon@phenoData$Group == "PND8"] <- "PND8"
salmon@phenoData$Batch <- "Lab"
salmon@phenoData$Batch[salmon@phenoData$Group %in% c("PND14", "PNW8")] <- "Lit"


# Save salmon data as R object ----
save(salmon, file = "./output/SC_controls_rnaseq_lab_lit_salmon.tds.RData", compress = T, compression_level = 3)


# SessionInfo ----
cat("# SessionInfo\n")
devtools::session_info()