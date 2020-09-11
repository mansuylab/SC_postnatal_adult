salmonImporter1 <- function(folder = NULL, files = NULL, anno = NULL, typeColumn = "type2",
                           geneColumn = "symbol", aggregationFilter = NULL) {
  if ((is.null(folder) & is.null(files)) | (!is.null(folder) &
    !is.null(files))) {
    stop("Exactly one of `folder` or `files` must be given.")
  }
  if (!is.null(folder)) {
    files <- list.files(folder,
      pattern = "quant.sf", recursive = T,
      full.names = T
    )
    message(paste("Qt file found:", files, collapse = "\n"))
  }
  if (is.null(anno)) {
    x <- fread(files[[1]], colClasses = c(
      "character", "integer",
      "numeric", "numeric", "numeric"
    ))$Name
    if (length(strsplit(x[[1]], "|", fixed = T)[[1]]) >=
      8) {
      anno <- t(sapply(x, FUN = function(x) {
        x <- strsplit(x, "|", fixed = T)[[1]]
        if (length(x) == 1) {
          x <- c(rep(x, 3), NA)
        }
        else {
          x <- x[c(1, 2, 6, 8)]
        }
      }))
      row.names(anno) <- anno[, 1]
      anno <- as.data.frame(anno[, 2:4], stringAsFactors = F)
      colnames(anno) <- c("gid", "symbol", "type")
    }
    else {
      stop("`anno` was not given, and the identifiers do not appear to be gencode records.")
    }
  }
  else {
    if (is.character(anno) && grepl("rda$|rdata$", anno,
      ignore.case = T
    )) {
      anno <- get(load(anno))
    }
  }
  snames <- basename(gsub("quant.sf$", "", files))
  library(data.table)
  library(jsonlite)
  cmdfiles <- lapply(gsub("quant.sf", "cmd_info.json", files),
    FUN = function(x) {
      try(lapply(read_json(x), toString), silent = T)
    }
  )
  hasLog <- !any(sapply(cmdfiles, class = "try-error", FUN = is))
  if (!hasLog) {
    warning("Could not find salmon log files for all samples; will skip importing method info.")
  }
  else {
    cmd <- rbindlist(cmdfiles, fill = TRUE)
    row.names(cmd) <- snames
    if (length(unique(cmd$salmon_version)) > 1) {
      warning("The different samples were quantified using different salmon versions! Will attempt aggregating them nonetheless, but the quantification might not be comparable.")
    }
    if (length(unique(cmd$index)) > 1) {
      warning("The different samples were quantified using different indexes! Will attempt aggregating them nonetheless, but the quantification might not be comparable.")
    }
    if (length(unique(cmd$libType)) > 1) {
      warning("The different samples were quantified using different libType parameters.")
    }
    libst <- lapply(gsub("quant.sf", "cmd_info.json", files),
      FUN = function(x) {
        lapply(read_json(x), toString)
      }
    )
    libst <- rbindlist(lapply(libst, FUN = function(x) x[-1]), fill = TRUE)
    row.names(cmd) <- snames
    a <- as.character(sapply(gsub("quant.sf", "logs/salmon_quant.log",
      files,
      fixed = T
    ), FUN = function(x) {
      x <- system(paste("grep \"Mapping rate\"", x), intern = T)
      x[length(x)]
    }))
    libst$output <- NULL
    libst$auxDir <- NULL
    libst$mapping_rate <- as.numeric(gsub("%", "", sapply(strsplit(a,
      " ",
      fixed = T
    ), FUN = function(x) {
      x[length(x)]
    }), fixed = T))
  }
  q <- lapply(files, colClasses = c(
    "character", "integer",
    "numeric", "numeric", "numeric"
  ), FUN = fread)
  names(q) <- snames
  if (length(unique(apply(sapply(q, dim), 2,
    collapse = " ",
    paste
  ))) > 1) {
    stop("The different quantification files did not use the same annotation and cannot be aggregated.")
  }
  els <- sapply(q, FUN = function(x) {
    x$EffectiveLength
  })
  ti <- data.frame(
    row.names = q[[1]]$Name, length = q[[1]]$Length,
    medianEffectiveLength = round(apply(els, 1, FUN = median)),
    effectiveLengthSD = round(apply(els, 1, FUN = sd))
  )
  ti <- cbind(ti, anno[row.names(ti), ])
  counts <- round(sapply(q, FUN = function(x) {
    x$NumReads
  }), 1)
  tpm <- round(sapply(q, FUN = function(x) {
    x$TPM
  }), 2)
  row.names(counts) <- row.names(ti)
  row.names(tpm) <- row.names(ti)
  rm(q)
  cc <- as.list(match.call())
  for (f in c("typeColumn", "geneColumn", "aggregationFilter")) {
    if (!(f %in% names(cc))) {
      cc[[f]] <- get(f)
    }
  }
  o <- new("transcriptomicDataset",
    call = as.call(cc), cmd = data.frame(row.names = colnames(counts)),
    libStats = data.frame(row.names = colnames(counts)),
    phenoData = data.frame(row.names = colnames(counts)),
    tx.info = ti, tx.counts = as.matrix(counts), tx.tpm = as.matrix(tpm)
  )
  message("Finished importing the data, now aggregation at the biotype and gene level.")
  o <- plgINS:::.transcriptomicDataset.aggregation(o,
    typeColumn = typeColumn,
    geneColumn = geneColumn, aggregationFilter = aggregationFilter
  )
  o@norm.factors <- calcNormFactors(o@gene.counts)
  if (hasLog) {
    o@libStats <- as.data.frame(libst)
    o@cmd <- as.data.frame(cmd)
  }
  return(o)
}
