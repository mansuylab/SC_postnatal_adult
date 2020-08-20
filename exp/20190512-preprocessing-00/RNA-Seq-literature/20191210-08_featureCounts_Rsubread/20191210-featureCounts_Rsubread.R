#!/usr/bin/env Rscript

# Library ----
library(Rsubread)

# allow arguments from command line
args <- commandArgs(trailingOnly = T)

# Parameters here were chosen to be in harmony with the options suggested by PLG in Issue #4 of SC_longRNA:
# primaryOnly=TRUE corresponds to "--primary"
# strandSpecific = 2 corresponds to "-s 2"
# useMetaFeatures=FALSE corresponds to "-f"
# allowMultiOverlap=TRUE corresponds to "-O"
fc <- featureCounts(
  files = args[1], primaryOnly = TRUE, strandSpecific = 2, useMetaFeatures = FALSE,
  allowMultiOverlap = TRUE, nthreads = args[2],
  annot.inbuilt = "mm10", annot.ext = "./input/gencode.vM18.chr_patch_hapl_scaff.annotation.gtf.gz",
  countMultiMappingReads = FALSE,
  isGTFAnnotationFile = T, GTF.featureType = "exon", GTF.attrType = "gene_id"
)

# save the featureCounts matrix for later use in DEXSeq script
name <- gsub(pattern = ".bam", replacement = "", x = basename(args[1]))
save(fc, file = paste0("./output/", name, ".RData"), compress = T, compression_level = 3)

cat("# SessionInfo\n\n")
devtools::session_info()