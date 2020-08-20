.write_bed <- function (x, path) 
{
  readr::write_tsv(x, path, col_names = FALSE)
}

findMotifs <- function(x, path, genome, motif_length = c(8, 10, 12), scan_size = 100,
                       optimize_count = 8, background = "automatic", local_background = FALSE,
                       only_known = FALSE, only_denovo = FALSE, fdr_num = 0, 
                       motif_file,
                       cores = parallel::detectCores(),
                       cache = .calc_free_mem() / 4, overwrite = FALSE, keep_minimal = FALSE) {
  if (overwrite == FALSE & dir.exists(path)) {
    stop("Output directory exists (set `overwrite = TRUE` to bypass)")
  }
  if (background != "automatic" && local_background != FALSE) {
    stop("`background` and `local_background` are mutually exclusive; use only one")
  }
  if (only_known != FALSE & only_denovo != FALSE) {
    stop("Both `only_known` and `only_denovo` set to `TRUE`; pick one")
  }
  if ("data.frame" %in% class(x)) {
    target_bed <- tempfile("target_")
    .write_bed(x, path = target_bed)
  }
  else {
    if (file.exists(x) != TRUE) {
      stop("Check that your bed file for `x` exists")
    }
    target_bed <- x
  }
  if (!("automatic" %in% background)) {
    if ("data.frame" %in% class(background)) {
      background_bed <- tempfile("background_")
      .write_bed(background, path = background_bed)
    }
    else {
      if (file.exists(background) != TRUE) {
        stop("Check that your bed file for `background` exists")
      }
      background_bed <- background
    }
  }
  system(paste("mkdir -p", path))
  homer_base <- get_homer_bin()
  cmd <- paste(
    paste0(homer_base, "findMotifsGenome.pl"), target_bed,
    genome, path, "-len", paste0(motif_length, collapse = ","),
    "-size", scan_size, "-S", optimize_count, "-p", cores,
    "-cache", cache, "-fdr", fdr_num
  )
  if (!("automatic" %in% background)) {
    cmd <- paste(cmd, "-bg", background_bed)
  }
  if (local_background != FALSE) {
    cmd <- paste(cmd, "-local", local_background)
  }
  if (only_known == TRUE) {
    cmd <- paste(cmd, "-nomotif")
  }
  if (only_denovo == TRUE) {
    cmd <- paste(cmd, "-noknown")
  }
  if (scan_size == "given") {
    cmd <- paste(cmd, "-chopify")
  }
  if (!is.null(motif_file)) {
    cmd <- paste(cmd, "-mknown", motif_file)
  }
  system(cmd)
  if (keep_minimal == TRUE) {
    extra_files <- c(
      "homerResults.html", "knownResults.html",
      "homerMotifs.motifs*", "motifFindingParameters.txt",
      "seq.autonorm.tsv", "*tmp*"
    )
    extra_dirs <- c("homerResults", "knownResults", "randomizations")
    remove_extra <- paste(c(
      paste0("rm -f ", path, "/", extra_files),
      paste0("rm -Rf ", path, "/", extra_dirs)
    ), collapse = "; ")
    system(remove_extra)
  }
  system("rm -f *.tmp")
}
