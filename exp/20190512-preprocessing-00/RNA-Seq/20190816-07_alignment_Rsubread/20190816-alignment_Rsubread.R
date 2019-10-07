#!/usr/bin/env Rscript

# load library
library(Rsubread)

# allow arguments from command line
args <- commandArgs(trailingOnly = T)

out <- gsub(pattern = "_trimmed.fq.gz", replacement = "", x = basename(args[1]))
outFile <- paste0("./output/", out, ".bam")

align.stat <- subjunc(
  index = "./input/Rsubread/GRCm38", readfile1 = args[1], readfile2 = NULL,
  input_format = "gzFASTQ", output_format = "BAM",
  output_file = outFile, phredOffset = 33,
  maxMismatches = 2, unique = FALSE, nBestLocations = 1,
  indels = 2, PE_orientation = "fr", nthreads = parallel::detectCores() - 1,
  sortReadsByCoordinates = T,
  useAnnotation = T, annot.inbuilt = "mm10",
  annot.ext = "./input/gencode.vM18.chr_patch_hapl_scaff.annotation.gtf.gz",
  isGTF = T, GTF.featureType = "exon", GTF.attrType = "gene_id"
)

cat(paste0(propmapped(outFile), "\n"))

cat(paste0(devtools::session_info(), "\n"))
