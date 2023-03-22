library(ComplexHeatmap)

### Set variables
outdir <- "~"

source(file.path(outdir, "Figure 3 - Fragment Ratio.R"))
source(file.path(outdir, "Figure 3 - Fragment Ratio PRC1.R"))

Figure <- Heatmap %v% PRC1_heatmap

pdf(file.path(outdir, "Figure 3.pdf"), width = 10, height = 9.5)
draw(Figure, merge_legend = TRUE, ht_gap = unit(0.25, "cm"))
dev.off()
