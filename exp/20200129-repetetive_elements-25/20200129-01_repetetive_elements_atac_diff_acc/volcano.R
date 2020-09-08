# Library ----
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
ht_opt$message <- FALSE



# Colors ----
cols <- magma(n = 10)[c(3, 6, 9)]
cols <- c(PND8 = cols[3], PND15 = cols[2], Adult = cols[1])


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
    size = 10,
    face = "plain",
    vjust = 1
  ),
  
  # axis text
  axis.text.x = element_text(
    angle = 0,
    size = 6,
    # vjust = 0.5
  ),
  axis.text.y = element_text(
    angle = 0,
    size = 6,
    # vjust = 0.5
  ),
  axis.title = element_text(
    size = 8,
    face = "plain",
  ),
  
  # legend
  legend.position = "top",
  legend.key = element_blank(),
  legend.key.size = unit(0, "cm"),
  legend.text = element_text(
    size = 8,
    margin = margin(r = 0, t = 0, b = 0, l = 0, unit = "pt"),
    vjust = 0.5, hjust = 0.5
  ),
  legend.title = element_text(size = 8, face = "bold", margin = margin(l = 0, unit = "pt")),
  legend.justification = "center",
  legend.margin = margin(c(0, 0, 0, 0)),
  legend.box.margin = margin(0, 0, -7, 0),
  legend.key.height = unit(0, "cm"),
  legend.key.width = unit(0, "cm"),
  legend.spacing = unit(1, "mm"),
  plot.margin = grid::unit(c(1, 1, 1, 1), "mm"),
  strip.background = element_rect(fill = NA, colour = "black"),
  strip.text.x = element_text(face = "bold", size = 8)
)



# Volcano plot ----
tab$Significant <- "Not Significant"
tab$Significant[tab$qvalue <= 0.05 & abs(tab$logFC) >= 1] <- "Negative"

tab$Significant[tab$Significant == "Negative"] <- ifelse(tab$logFC[tab$Significant == "Negative"] > 0, "Positive", "Negative")

tab$Significant[tab$Significant == "Negative"] <- "Less acc. in Adult"
tab$Significant[tab$Significant == "Positive"] <- "More acc. in Adult"

tab$Significant <- factor(x = tab$Significant, levels = unique(tab$Significant))

an <- split(x = tab, f = tab$class_id)
an1 <- lapply(an, function(x) {
  return(data.frame(
    x = c(-5, 5),
    y = c(5, 5),
    label = c(
      paste(
        "n =",
        nrow(x[x$logFC <= -1 & x$qvalue <= 0.05, ])
      ),
      paste(
        "n =",
        nrow(x[x$logFC >= 1 & x$qvalue <= 0.05, ])
      )
    ),
    col = c("#BC3C29FF", "#0072B5FF")
  ))
})

anno <- ldply(an1, data.frame, .id = "class_id")
anno <- subset(anno, class_id %in% c("LTR", "LINE", "SINE"))


a <- ggplot(subset(data.frame(tab), class_id %in% c("LTR", "LINE", "SINE")), aes(x = logFC, y = -log10(qvalue))) +
  geom_point(aes(color = Significant), size = 2, alpha = 0.5) +
  scale_color_manual(values = c("black", "#BC3C29FF", "#0072B5FF")) +
  xlim(
    -(ceiling(max(abs(tab$logFC)))),
    ceiling(max(abs(tab$logFC)))
  ) +
  ylim(0, max(-log10(tab$qvalue), na.rm = TRUE) + 1) +
  xlab(bquote(~ Log[2] ~ "fold change (FC)")) +
  ylab(bquote(~ -Log[10] ~ adjusted ~ italic(P))) +
  labs(
    title = "",
    subtitle = "", caption = paste0("Total = ", scales::comma(nrow(tab)), " tested TEs")
  ) + theme_bw(base_family = "Helvetica") +
  th +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "longdash",
    colour = "black",
    size = 0.4
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "longdash",
    colour = "black",
    size = 0.4
  ) +
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(
        size = 4
      )
    )
  ) +
  scale_shape_manual(
    values = c(
      `Not Significant` = 19,
      Positive = 19,
      Negative = 19
    ),
    labels = c(
      `Not Significant` = "Not significant",
      Positive = "More acc. in Adult",
      Negative = "Less acc. in Adult"
    )
  ) +
  facet_wrap(~class_id)

ag <- egg::set_panel_size(p = a, width = unit(5.5, "cm"), height = unit(6, "cm"))
ggsave(plot = ag, filename = "output/a.png", width = unit(8.5, "in"), height = unit(11, "in"))
