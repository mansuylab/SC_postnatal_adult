library(markdown)
library(knitr)
library(EnrichedHeatmap)
library(GetoptLong)
library(circlize)
library(RColorBrewer)

download.file("https://jokergoo.github.io/supplementary/EnrichedHeatmap-supplementary/roadmap_normalized_matrices.RData",
  destfile = "roadmap_normalized_matrices.RData"
)
load("roadmap_normalized_matrices.RData")
file.remove("roadmap_normalized_matrices.RData")

SAMPLE
COLOR

## ---- eval = FALSE--------------------------------------------------------------------------------
#  # this chunk of code is only for demonstration
#  mat_corr = normalizeToMatrix(cr, tss, mapping_column = "gene_id", value_column = "corr",
#      mean_mode = "absolute", ...)

#  # this chunk of code is only for demonstration
#  mat_neg_cr = normalizeToMatrix(sig_neg_cr, tss, mapping_column = "gene_id",
#      mean_mode = "absolute", ...)
#  # this chunk of code is only for demonstration
#  meth_mean = meth
#  mcols(meth_mean) = data.frame(mean_meth = rowMeans(mcols(meth)))
#  meth_mat_mean = normalizeToMatrix(meth_mean, tss, value_column = "mean_meth",
#      mode = "absolute", ...)

#  # this chunk of code is only for demonstration
#  meth_mean_1 = meth
#  mcols(meth_mean_1) = data.frame(mean_meth = rowMeans(mcols(meth[, SAMPLE$subgroup == "subgroup1"])))
#  meth_mat_mean_1 = normalizeToMatrix(meth_mean_1, tss, value_column = "mean_meth",
#      mode = "absolute", ...)
#  meth_mean_2 = meth
#  mcols(meth_mean_2) = data.frame(mean_meth = rowMeans(mcols(meth[, SAMPLE$subgroup == "subgroup2"])))
#  meth_mat_mean_2 = normalizeToMatrix(meth_mean_2, tss, value_column = "mean_meth",
#      mode = "absolute", ...)
#  meth_mat_diff = meth_mat_mean_1 - meth_mat_mean_2

## -------------------------------------------------------------------------------------------------
length(tss)

## ---- eval = FALSE--------------------------------------------------------------------------------
#  # this chunk of code is only for demonstration
#  for(i in n_sample) {
#      # assume the column name for the signals is called 'density'
#      hm_list[[i]] = normalizeToMatrix(peak[[i]], tss, value_column = "density", keep = c(0, 0.99))
#  }

## ---- eval = FALSE--------------------------------------------------------------------------------
#  # this chunk of code is only for demonstration
#  hm_mat_mean = getSignalsFromList(hm_list, mean)

## ---- eval = FALSE--------------------------------------------------------------------------------
#  # this chunk of code is only for demonstration
#  # hm_list_1 and hm_list_2 are normalized matrices for subgroup1 and subgroup2 separatedly
#  hm_mat_mean_1 = getSignalsFromList(hm_list_1, mean)
#  hm_mat_mean_2 = getSignalsFromList(hm_list_2, mean)
#  hm_mat_diff = hm_mat_mean_1 - hm_mat_mean_2

## ---- eval = FALSE--------------------------------------------------------------------------------
#  # this chunk of code is only for demonstration
#  hm_corr = getSignalsFromList(hm_list, fun = function(x, i) {
#      cor(x, expr[i, ], method = "spearman") # x = array[i, j, ]
#  })

## ---- eval = FALSE--------------------------------------------------------------------------------
#  # this chunk of code is only for demonstration
#  mat_cgi = normalizeToMatrix(cgi, tss, mean_mode = "absolute", ...)

expr_mean <- rowMeans(expr[, SAMPLE$subgroup == "subgroup1"]) -
  rowMeans(expr[, SAMPLE$subgroup == "subgroup2"])
expr_split <- ifelse(expr_mean > 0, "high", "low")
expr_split <- factor(expr_split, levels = c("high", "low"))

set.seed(123)
upstream_index <- length(attr(meth_mat_mean, "upstream_index"))
meth_split <- kmeans(meth_mat_mean[, seq(round(upstream_index * 0.8), round(upstream_index * 1.4))],
  centers = 2
)$cluster
x <- tapply(
  rowMeans(meth_mat_mean[, seq(round(upstream_index * 0.8), round(upstream_index * 1.4))]),
  meth_split, mean
)
od <- structure(order(x), names = names(x))
meth_split <- paste0("cluster", od[as.character(meth_split)])

combined_split <- paste(meth_split, expr_split, sep = "|")

tb <- table(combined_split)
tb
tb["cluster2|high"] / sum(tb)

l <- combined_split != "cluster2|high"
tss <- tss[l]
expr <- expr[l, ]
hist_mat_corr_list <- lapply(hist_mat_corr_list, function(x) x[l, ])
hist_mat_mean_list <- lapply(hist_mat_mean_list, function(x) x[l, ])
hist_mat_diff_list <- lapply(hist_mat_diff_list, function(x) x[l, ])
mat_neg_cr <- mat_neg_cr[l, ]
mat_cgi <- mat_cgi[l, ]
meth_mat_corr <- meth_mat_corr[l, ]
meth_mat_mean <- meth_mat_mean[l, ]
meth_mat_diff <- meth_mat_diff[l, ]
expr_split <- expr_split[l]
meth_split <- meth_split[l]
combined_split <- combined_split[l]
n_row_cluster <- length(unique(combined_split))

merge_row_order <- function(l_list) {
  do.call("c", lapply(l_list, function(l) {
    if (sum(l) == 0) {
      return(integer(0))
    }
    if (sum(l) == 1) {
      return(which(l))
    }
    dend1 <- as.dendrogram(hclust(dist_by_closeness(mat_neg_cr[l, ])))
    dend1 <- reorder(dend1, -enriched_score(mat_neg_cr[l, ]))
    od <- order.dendrogram(dend1)
    which(l)[od]
  }))
}

row_order <- merge_row_order(list(
  combined_split == "cluster1|high",
  combined_split == "cluster1|low",
  combined_split == "cluster2|low"
))

dend1 <- as.dendrogram(hclust(dist(t(expr[, SAMPLE$subgroup == "subgroup1"]))))
hc1 <- as.hclust(reorder(dend1, colMeans(expr[, SAMPLE$subgroup == "subgroup1"])))
expr_col_od1 <- hc1$order
dend2 <- as.dendrogram(hclust(dist(t(expr[, SAMPLE$subgroup == "subgroup2"]))))
hc2 <- as.hclust(reorder(dend2, colMeans(expr[, SAMPLE$subgroup == "subgroup2"])))
expr_col_od2 <- hc2$order
expr_col_od <- c(
  which(SAMPLE$subgroup == "subgroup1")[expr_col_od1],
  which(SAMPLE$subgroup == "subgroup2")[expr_col_od2]
)

ht_list <- Heatmap(expr,
  name = "expr", show_row_names = FALSE,
  show_column_names = FALSE, width = unit(4, "cm"), show_column_dend = FALSE,
  cluster_columns = FALSE, column_order = expr_col_od,
  top_annotation = HeatmapAnnotation(
    df = SAMPLE[, -1], col = COLOR,
    annotation_name_side = "left"
  ),
  column_title = "Expression", column_title_gp = gpar(fontsize = 12),
  show_row_dend = FALSE, use_raster = TRUE
)

library(genefilter)
df <- rowttests(expr, factor(SAMPLE$subgroup))
top_genes <- rownames(df[order(df$p.value)[1:20], ])

index <- which(rownames(expr) %in% top_genes)
labels <- gene_symbol[rownames(expr)[index]]
ht_list <- rowAnnotation(sig_gene = anno_mark(
  at = index, labels = labels,
  side = "left", labels_gp = gpar(fontsize = 10),
  extend = unit(c(1, 0), "cm")
)) + ht_list

gl <- width(gene[names(tss)])
gl[gl > quantile(gl, 0.95)] <- quantile(gl, 0.95)
ht_list <- ht_list + rowAnnotation(gene_len = anno_points(gl,
  size = unit(1, "mm"),
  gp = gpar(col = "#00000040"),
  axis_param = list(at = c(0, 1e5, 2e5), labels = c("0bp", "100bp", "200bp")),
  width = unit(1.5, "cm")
))

axis_name <- c("-5kb", "TSS", "10kb")
ht_list <- ht_list + EnrichedHeatmap(mat_cgi,
  col = c("white", "darkorange"), name = "CGI",
  column_title = "CGI", column_title_gp = gpar(fontsize = 12),
  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(
    col = "darkorange",
    lty = 1:n_row_cluster
  ), axis_param = list(side = "right", facing = "inside"))),
  axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), use_raster = TRUE
)

bg_col <- brewer.pal(8, "Set2")
cor_col_fun <- colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))
ht_list <- ht_list + EnrichedHeatmap(meth_mat_corr,
  col = cor_col_fun, name = "meth_corr",
  top_annotation = HeatmapAnnotation(lines = anno_enriched(
    gp = gpar(
      pos_col = "red",
      neg_col = "darkgreen", lty = 1:n_row_cluster
    ),
    axis_param = list(side = "right", facing = "inside")
  )),
  column_title = "meth_corr", column_title_gp = gpar(fontsize = 12, fill = bg_col[1]),
  axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), use_raster = TRUE
)

meth_col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
ht_list <- ht_list + EnrichedHeatmap(meth_mat_mean,
  col = meth_col_fun, name = "meth_mean",
  column_title = "meth_mean", column_title_gp = gpar(fontsize = 12, fill = bg_col[1]),
  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(
    col = "red",
    lty = 1:n_row_cluster
  ), axis_param = list(side = "right", facing = "inside"))),
  axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), use_raster = TRUE
)

generate_diff_color_fun <- function(x) {
  q <- quantile(x, c(0.05, 0.95))
  max_q <- max(abs(q))
  colorRamp2(c(-max_q, 0, max_q), c("#3794bf", "#FFFFFF", "#df8640"))
}

ht_list <- ht_list + EnrichedHeatmap(meth_mat_diff,
  name = "meth_diff",
  col = generate_diff_color_fun(meth_mat_diff),
  column_title = "meth_diff", column_title_gp = gpar(fontsize = 12, fill = bg_col[1]),
  top_annotation = HeatmapAnnotation(lines = anno_enriched(
    gp = gpar(
      pos_col = "#df8640",
      neg_col = "#3794bf", lty = 1:n_row_cluster
    ),
    axis_param = list(side = "right", facing = "inside")
  )),
  axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), use_raster = TRUE
)

ht_list_2 <- NULL
ht_list_1 <- NULL
mark_name <- names(hist_mat_corr_list)
for (i in seq_along(hist_mat_corr_list)) {
  # heatmaps for the 2nd, 3th and 4th histone modifications are assigned to a new `ht_list`
  if (i == 2) {
    ht_list_1 <- ht_list
    ht_list <- NULL
  }

  ht_list <- ht_list + EnrichedHeatmap(hist_mat_corr_list[[i]],
    col = cor_col_fun,
    name = qq("@{mark_name[i]}_corr"), column_title = qq("@{mark_name[i]}_corr"),
    column_title_gp = gpar(fontsize = 12, fill = bg_col[i + 1]),
    top_annotation = HeatmapAnnotation(lines = anno_enriched(
      gp = gpar(
        pos_col = "red",
        neg_col = "darkgreen", lty = 1:n_row_cluster
      ),
      axis_param = list(side = "right", facing = "inside")
    )),
    axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), use_raster = TRUE
  )

  ht_list <- ht_list + EnrichedHeatmap(hist_mat_mean_list[[i]],
    col = colorRamp2(c(0, quantile(hist_mat_mean_list[[i]], 0.95)), c("white", "purple")),
    name = qq("@{mark_name[i]}_mean"), column_title = qq("@{mark_name[i]}_mean"),
    column_title_gp = gpar(fontsize = 12, fill = bg_col[i + 1]),
    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(
      col = "purple",
      lty = 1:n_row_cluster
    ), axis_param = list(side = "right", facing = "inside"))),
    axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), use_raster = TRUE
  )

  ht_list <- ht_list + EnrichedHeatmap(hist_mat_diff_list[[i]],
    col = generate_diff_color_fun(hist_mat_diff_list[[i]]),
    name = qq("@{mark_name[i]}_diff"), column_title = qq("@{mark_name[i]}_diff"),
    column_title_gp = gpar(fontsize = 12, fill = bg_col[i + 1]),
    top_annotation = HeatmapAnnotation(lines = anno_enriched(
      gp = gpar(
        pos_col = "#df8640",
        neg_col = "#3794bf", lty = 1:n_row_cluster
      ),
      axis_param = list(side = "right", facing = "inside")
    )),
    axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), use_raster = TRUE
  )
}
ht_list_2 <- ht_list

split <- as.vector(combined_split)
split[combined_split == "cluster1|high"] <- "cluster1"
split[combined_split == "cluster1|low"] <- "cluster2"
split[combined_split == "cluster2|low"] <- "cluster3"

ht_list_1 <- Heatmap(expr_split,
  show_row_names = FALSE, name = "expr_diff",
  col = c("high" = "red", "low" = "darkgreen"),
  show_column_names = FALSE, width = unit(2, "mm")
) + ht_list_1

ht_list_1 <- draw(ht_list_1,
  cluster_rows = FALSE, row_order = row_order, show_row_dend = FALSE,
  row_split = split, heatmap_legend_side = "bottom", ht_gap = unit(2, "mm")
)

add_boxplot_of_gene_length <- function(ht_list) {
  row_order_list <- row_order(ht_list)
  lt <- lapply(row_order_list, function(ind) gl[ind])
  bx <- boxplot(lt, plot = FALSE)$stats
  n <- length(row_order_list)

  decorate_annotation("gene_len", slice = 1, {
    rg <- range(bx)
    rg[1] <- rg[1] - (rg[2] - rg[1]) * 0.1
    rg[2] <- rg[2] + (rg[2] - rg[1]) * 0.1
    pushViewport(viewport(
      y = unit(1, "npc") + unit(1, "mm"), just = "bottom",
      height = unit(2, "cm"), yscale = rg, xscale = c(0.5, n + 0.5)
    ))
    grid.rect()
    for (i in 1:n) {
      grid.boxplot(pos = i, lt[[i]], gp = gpar(lty = i), outline = FALSE)
    }
    grid.text("Gene length",
      y = unit(1, "npc") + unit(2.5, "mm"),
      gp = gpar(fontsize = 12), just = "bottom"
    )
    upViewport()
  })
}

add_boxplot_of_gene_length(ht_list_1)


ht_list_2 <- Heatmap(expr_split,
  show_row_names = FALSE, name = "expr_diff",
  col = c("high" = "red", "low" = "darkgreen"),
  show_column_names = FALSE, width = unit(2, "mm")
) + ht_list_2
ht_list_2 <- draw(ht_list_2,
  cluster_rows = FALSE, row_order = row_order, show_row_dend = FALSE,
  row_split = split, heatmap_legend_side = "bottom", ht_gap = unit(2, "mm")
)


add_anno_enriched <- function(ht_list, name, ri, ci) {
  pushViewport(viewport(layout.pos.row = ri, layout.pos.col = ci))
  extract_anno_enriched(ht_list, name, newpage = FALSE)
  upViewport()
}

pushViewport(viewport(layout = grid.layout(nr = 3, nc = 6)))
add_anno_enriched(ht_list_1, "meth_corr", 1, 1)
add_anno_enriched(ht_list_1, "meth_mean", 1, 2)
add_anno_enriched(ht_list_1, "meth_diff", 1, 3)
add_anno_enriched(ht_list_1, "CGI", 1, 4)
add_anno_enriched(ht_list_1, "H3K4me3_corr", 2, 1)
add_anno_enriched(ht_list_1, "H3K4me3_mean", 2, 2)
add_anno_enriched(ht_list_1, "H3K4me3_diff", 2, 3)
add_anno_enriched(ht_list_2, "H3K4me1_corr", 2, 4)
add_anno_enriched(ht_list_2, "H3K4me1_mean", 2, 5)
add_anno_enriched(ht_list_2, "H3K4me1_diff", 2, 6)
add_anno_enriched(ht_list_2, "H3K27ac_corr", 3, 1)
add_anno_enriched(ht_list_2, "H3K27ac_mean", 3, 2)
add_anno_enriched(ht_list_2, "H3K27ac_diff", 3, 3)
add_anno_enriched(ht_list_2, "H3K27me3_corr", 3, 4)
add_anno_enriched(ht_list_2, "H3K27me3_mean", 3, 5)
add_anno_enriched(ht_list_2, "H3K27me3_diff", 3, 6)
upViewport()