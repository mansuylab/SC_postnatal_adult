#' annoPeaks
#'
#' Annotates a GRanges on the basis of a (GENCODE-like) gtf file.
#' @author Pierre-Luc Germain germain@hifo.uzh.ch
#' @param peaks A GRanges object
#' @param gtf A GRanges object of a GENCODE-like gtf, or the path to such a gtf
#'  file, or a TxDb object.
#' @param proximal The threshold(s) for TSS proximal regions
#'
#' @return The sorted `peaks` object with additional annotation columns.
#' @export
annoPeaks <- function(peaks, gtf, proximal = c(2500, 1000)) {
  library(GenomicRanges)
  peaks <- sort(peaks)
  if (is.character(gtf) && length(gtf) == 1) {
    gtf <- rtracklayer::import.gff(gtf)
  }
  if (is(gtf, "TxDb") || is(gtf, "EnsDb")) {
    if (is(gtf, "TxDb")) {
      tx <- transcripts(gtf, columns = c("tx_name", "gene_id"))
      tx$transcript_id <- tx$tx_name
      tx$tx_name <- NULL
      tx$gene_name <- tx$gene_id
    } else {
      tx <- transcripts(gtf, columns = c("tx_id", "gene_id", "gene_name"))
      tx$transcript_id <- tx$tx_id
      tx$tx_id <- NULL
    }
    ex <- exons(gtf, columns = NULL)
    tx$type <- "transcript"
    mcols(ex)$type <- "exon"
    gtf <- sort(c(tx, ex))
    gtf$type <- as.factor(gtf$type)
  }
  if (!is(gtf, "GRanges") || !is(peaks, "GRanges")) {
    stop("`gtf` and `peaks` should be GRanges objects.")
  }
  seqlevelsStyle(gtf) <- seqlevelsStyle(peaks)
  gtf <- gtf[gtf$type %in% c("transcript", "exon", "dispersed_repeat")]
  sl <- intersect(seqlevels(gtf), seqlevels(peaks))
  if (length(sl) == 0) stop("No seqlevel in common!")
  gtf <- keepSeqlevels(gtf, sl, pruning.mode = "coarse")
  peaks <- keepSeqlevels(peaks, sl, pruning.mode = "coarse")
  gtf$type <- relevel(droplevels(gtf$type), "exon")
  tss <- gtf[gtf$type == "transcript"]
  tss1 <- tss[strand(tss) == "+"]
  tss2 <- tss[strand(tss) == "-"]
  tss <- c(
    GRanges(seqnames(tss1), IRanges(start(tss1), width = 1),
      strand = strand(tss1),
      mcols(tss1)[, c("transcript_id", "gene_id", "gene_name")]
    ),
    GRanges(seqnames(tss2), IRanges(end(tss2), width = 1),
      strand = strand(tss2),
      mcols(tss2)[, c("transcript_id", "gene_id", "gene_name")]
    )
  )
  d <- distanceToNearest(peaks, tss)
  peaks$distance2nearestTSS <- NA_integer_
  peaks$distance2nearestTSS[d@from] <- mcols(d)$distance
  tmp <- cbind(start(peaks)[d@from], end(peaks)[d@from]) - start(tss)[d@to]
  peaks$distance2nearestTSS[d@from] <- apply(tmp, 1, FUN = function(x) {
    if (length(unique(sign(x))) > 1) {
      return(0)
    }
    x[which.min(abs(x))]
  })
  peaks$nearestTSS <- peaks$nearestTSS.gene_name <- NA_character_
  mcols(peaks)[d@from, c("nearestTSS", "nearestTSS.gene_name")] <- mcols(tss)[d@to, c("transcript_id", "gene_name")]
  o <- findOverlaps(peaks, gtf)
  ll <- sapply(split(o@to, o@from), FUN = function(x) sort(gtf$type[x])[1])
  peaks$overlap <- "intergenic"
  peaks$overlap[as.numeric(names(ll))] <- c("exonic", "intronic", "repeat")[as.numeric(ll)]
  peaks$class <- peaks$overlap
  proximal <- sort(proximal, decreasing = TRUE)
  for (i in seq_along(proximal)) {
    lab <- paste0(
      "proximal ",
      ifelse(i != length(proximal), paste0(">", proximal[i + 1], "&"), ""),
      "<=", proximal[i], "bp"
    )
    peaks$class[abs(peaks$distance2nearestTSS) <= proximal[i]] <- lab
  }
  peaks$class[abs(peaks$distance2nearestTSS) == 0] <- "TSS"
  peaks$overlap <- factor(peaks$overlap)
  peaks$class <- as.factor(peaks$class)
  peaks
}
