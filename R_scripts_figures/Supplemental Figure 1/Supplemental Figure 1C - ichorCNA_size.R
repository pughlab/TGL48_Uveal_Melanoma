library(tidyverse)
library(data.table)
library(ggpubr)

### Set paths
healthy <- "TGL48_0009_Ly_R_PE_677_WG.merged.sorted_deduped.correctedDepth.txt"
healthy_seg <-"TGL48_0009_Ly_R_PE_677_WG.merged.sorted_deduped.seg.txt"
ichor <- "TGL48_0009_Ct_T_PE_327_WG.merged.sorted_deduped.correctedDepth.txt"
ichor_seg <-"TGL48_0009_Ct_T_PE_327_WG.merged.sorted_deduped.seg.txt"
short <- "TGL48_0009_Ct_T_PE_327_WG.merged.sorted_deduped_short.correctedDepth.txt"
short_seg <-"TGL48_0009_Ct_T_PE_327_WG.merged.sorted_deduped_short.seg.txt"
outdir <- ""

### Read in data
healthy <- read.delim(healthy)
healthy_seg <- read.delim(healthy_seg)
ichor <- read.delim(ichor)
ichor_seg <- read.delim(ichor_seg)
short <- read.delim(short)
short_seg <- read.delim(short_seg)

### Format chromosome columns
armlevels <- c(1:22, "X")

healthy$chr <- factor(healthy$chr, levels = armlevels)
ichor$chr <- factor(ichor$chr, levels = armlevels)
short$chr <- factor(short$chr, levels = armlevels)

healthy_seg$chrom <- factor(healthy_seg$chrom, levels = armlevels)
ichor_seg$chrom <- factor(ichor_seg$chrom, levels = armlevels)
short_seg$chrom <- factor(short_seg$chrom, levels = armlevels)

### Merge seg info with ichorCNA copy number data
### Healthy
loc <- healthy[, 1:3]
seg_loc <- healthy_seg[, 2:4]
colnames(seg_loc) <- c("chr", "start", "end")

setDT(loc)
setDT(seg_loc)
setkey(seg_loc)
overlaps <- foverlaps(loc, seg_loc, type="within", nomatch=0L)
overlaps <- merge(overlaps, healthy_seg, by.x = c("chr", "start", "end"), by.y = c("chrom", "start", "end"))
healthy <- merge(healthy, overlaps, by.x = c("chr", "start", "end"), by.y = c("chr", "i.start", "i.end"))
rm(loc, seg_loc)

### All Fragments
loc <- ichor[, 1:3]
seg_loc <- ichor_seg[, 2:4]
colnames(seg_loc) <- c("chr", "start", "end")

setDT(loc)
setDT(seg_loc)
setkey(seg_loc)
overlaps <- foverlaps(loc, seg_loc, type="within", nomatch=0L)
overlaps <- merge(overlaps, ichor_seg, by.x = c("chr", "start", "end"), by.y = c("chrom", "start", "end"))
ichor <- merge(ichor, overlaps, by.x = c("chr", "start", "end"), by.y = c("chr", "i.start", "i.end"))
rm(loc, seg_loc)

### Short Fragments
loc <- short[, 1:3]
seg_loc <- short_seg[, 2:4]
colnames(seg_loc) <- c("chr", "start", "end")

setDT(loc)
setDT(seg_loc)
setkey(seg_loc)
overlaps <- foverlaps(loc, seg_loc, type="within", nomatch=0L)
overlaps <- merge(overlaps, short_seg, by.x = c("chr", "start", "end"), by.y = c("chrom", "start", "end"))
short <- merge(short, overlaps, by.x = c("chr", "start", "end"), by.y = c("chr", "i.start", "i.end"))
rm(loc, seg_loc)

### Set themes and plot layouts
mytheme <- theme_classic(base_size=12) + theme(
  axis.text.x = element_blank(),
  axis.ticks.x=element_blank(),
  strip.text.x = element_text(size=11),
  strip.text.y = element_text(size=12),
  axis.title.x = element_text(face="bold", size=17),
  axis.title.y = element_text(size=15),
  axis.text.y = element_text(size=15),
  plot.title = element_text(size=15),
  legend.position = "none",
  legend.title = element_text(size=10),
  legend.text = element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_rect(fill="white", color="white"),
  panel.spacing.x=unit(0.1, "lines"))

### Order ichorCNA data for plotting
healthy <- healthy[order(factor(healthy$chr, levels = armlevels),
                         healthy$start), ]
healthy$chr <- factor(healthy$chr, levels = armlevels)
healthy$bin <- c(1:nrow(healthy))

ichor <- ichor[order(factor(ichor$chr, levels = armlevels),
                         ichor$start), ]
ichor$chr <- factor(ichor$chr, levels = armlevels)
ichor$bin <- c(1:nrow(ichor))

short <- short[order(factor(short$chr, levels = armlevels),
                         short$start), ]
short$chr <- factor(short$chr, levels = armlevels)
short$bin <- c(1:nrow(short))

## Plot Fragmentation profile
g1 <- ggplot(healthy, aes(x=bin, y=log2_TNratio_corrected, color=Corrected_Call)) + 
  geom_point(size=0.75, alpha=0.75) + 
  geom_point(aes(x=bin, y=seg.median.logR, color=Corrected_Call), size = 1) +
  scale_color_manual(values = c("NEUT" = "#0000FF", "HETD" = "#006400", "GAIN" = "#8B0000", "AMP" = "#FF0000", "HLAMP" = "#FF0000")) +
  labs(x="", y="Lymphocytes", color="") + 
  facet_grid(~chr, switch="x",space="free_x", scales="free_x") + 
  coord_cartesian(xlim = NULL, ylim=c(-0.5,0.5), expand = TRUE) + 
  mytheme
g1

g2 <- ggplot(ichor, aes(x=bin, y=log2_TNratio_corrected, color=Corrected_Call)) + 
  geom_point(size=0.75, alpha=0.75) + 
  geom_point(aes(x=bin, y=seg.median.logR, color=Corrected_Call), size = 1) +
  scale_color_manual(values = c("NEUT" = "#0000FF", "HETD" = "#006400", "GAIN" = "#8B0000", "AMP" = "#FF0000", "HLAMP" = "#FF0000")) +
  labs(x="", y="All Fragments", color="") + 
  facet_grid(~chr, switch="x",space="free_x", scales="free_x") + 
  coord_cartesian(xlim = NULL, ylim=c(-0.5,0.5), expand = TRUE) + 
  mytheme
g2

g3 <- ggplot(short, aes(x=bin, y=log2_TNratio_corrected, color=Corrected_Call)) + 
  geom_point(size=0.75, alpha=0.75) + 
  geom_point(aes(x=bin, y=seg.median.logR, color=Corrected_Call), size = 1) +
  scale_color_manual(values = c("NEUT" = "#0000FF", "HETD" = "#006400", "GAIN" = "#8B0000", "AMP" = "#FF0000", "HLAMP" = "#FF0000")) +
  labs(x="", y="Short Fragments", color="") + 
  facet_grid(~chr, switch="x",space="free_x", scales="free_x") + 
  coord_cartesian(xlim = NULL, ylim=c(-0.5,0.5), expand = TRUE) + 
  mytheme
g3

p1 <- ggarrange(g1, g2, g3, ncol = 1, nrow = 3)
p1 <- annotate_figure(p1, top = text_grob("Patient 09 - 6 months", size = 18))
p1
ggsave(file.path(outdir, "Supplemental Figure 4 - ichorCNA_size.pdf"), p1, width=15, height=6, units="in")
