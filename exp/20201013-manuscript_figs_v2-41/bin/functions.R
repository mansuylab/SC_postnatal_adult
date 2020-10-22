# Libraries ----
library(tidyverse)
library(report)
library(ggpubr)
library(patchwork)
library(viridis)
library(ggfortify)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(plyr)
library(reshape2)
library(extrafont)
library(ggthemes)
library(EnhancedVolcano)
library(SummarizedExperiment)
library(csaw)
library(edgeR)
library(seriation)
library(ggsci)
library(RColorBrewer)
library(rcartocolor)
ht_opt$message <- FALSE


# Colors ----
cols <- rev(plasma(n = 10))[c(2,4,6,8,10)]
cols <- c(PND8 = cols[1], PND15 = cols[3], PND14 = cols[2], PNW8 = cols[4], Adult = cols[5])


# Colorbar ----
myPalette <- colorRampPalette(colors = viridis::plasma(n = 7)[c(2, 4, 6)])


# Theme ----
th <- theme(
  legend.background = element_rect(),
  
  # title, subtitle, and caption
  plot.title = element_text(
    angle = 0,
    size = 12,
    face = "bold",
    vjust = 1
  ),
  plot.subtitle = element_text(
    angle = 0,
    size = 10,
    face = "plain",
    vjust = 1
  ),
  plot.caption = element_text(
    angle = 0,
    size = 8,
    face = "plain",
    vjust = 1
  ),
  
  # axis text
  axis.text.x = element_text(
    angle = 0,
    size = 6,
    # vjust = 0.5,
    # hjust = 0.5
  ),
  axis.text.y = element_text(
    angle = 0,
    size = 6,
    # vjust = 0.5,
    # hjust = 0.5
  ),
  axis.title = element_text(
    size = 8,
    face = "plain",
  ),
  
  # legend
  legend.position = "top",
  legend.key = element_blank(),
  # legend.key.size = unit(0, "cm"),
  legend.text = element_text(
    size = 8,
    margin = margin(r = 0, t = 0, b = 0, l = 0, unit = "pt"),
    vjust = 0, hjust = 0
  ),
  legend.title = element_text(size = 8, face = "plain"),
  legend.justification = "bottom",
  legend.margin = margin(c(0, 0, 0, 0)),
  legend.box.margin = margin(0, 0, -7, 0),
  legend.key.height = unit(0, "cm"),
  legend.key.width = unit(0, "cm"),
  legend.spacing = unit(1, "mm"),
  plot.margin = grid::unit(c(0, 0, 0, 0), "mm"),
  strip.background = element_rect(fill = NA, colour = "black"),
  strip.text.x = element_text(face = "plain", size = 8)
)
