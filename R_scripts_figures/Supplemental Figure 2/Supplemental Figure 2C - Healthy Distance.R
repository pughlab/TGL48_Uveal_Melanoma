library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(matrixStats)

### Set variables
healthy_dir <- ""
outdir <- ""

data_ratio <- file.path(healthy_dir, "fragmentomics/HBC_fragment_ratios.txt")
data_samples <- file.path(healthy_dir, "HBC_sample_list.txt")
data_coverage <- file.path(healthy_dir, "coverage/HBC_coverage.txt")

### Read in data
data_ratio <- read.delim(data_ratio, check.names = FALSE)
data_samples <- read.delim(data_samples)
data_coverage <- read.delim(data_coverage)

### Set chromosomes
data_chr <- data_ratio[, c("seqnames", "arm", "start", "end")]
chr_names <- with(data_chr, paste0(seqnames, "_", start))
row.names(data_ratio) <- chr_names

### Format sample names
data_samples <- merge(data_samples, data_coverage, by = "sample")
data_samples$cohort <- factor(data_samples$cohort, levels = c("NF1_NCI", "Landau", "HCC", "Phallen"),
                              labels = c("Szymanski et al", "Zviran et al", "Wong et al", "Cristiano et al"))
data_samples <- data_samples[order(data_samples$cohort,
                                   data_samples$coverage), ]

### Calculate the healthy median and SD
data_healthy <- as.matrix(data_ratio[ , -c(1:4)])
data_healthy <- data_healthy[, data_samples$sample]
healthy_median <- rowMedians(data_healthy)
healthy_sd <- rowSds(data_healthy)

### Convert ratios to distance from healthy median
data_dist <- (data_healthy - healthy_median)/healthy_sd
data_dist <- as.matrix(t(data_dist))
data_dist <- data_dist[data_samples$sample, ]

### Calculate z-scores
normal_score <- rowSums(abs(data_dist), na.rm = TRUE)
normal_median <- median(normal_score)
normal_MAD <- mad(normal_score)
data_zscore <- abs((normal_score - normal_median)/normal_MAD)

### Set chr order
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
data_chr$arm <- factor(data_chr$arm, levels=armlevels)
data_chr$arm <- factor(data_chr$arm, levels=armlevels,
                       labels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                                  "7p","7q","8p","8q", "9p", "9q","","10q"," ","11q","  ",
                                  "12q","13q","14q","15q","   ","16q", "    ","17q","     ","18q",
                                  "      ", "       ", "        ", "         ","          ","           "))

### Set cohort annotation
data_cohort <- as.matrix(data_samples$cohort)
row.names(data_cohort) <- data_samples$sample

### Set coverage annotation
data_coverage <- as.matrix(data_samples$coverage)
row.names(data_coverage) <- data_samples$sample

## Set colours
col <- colorRamp2(c(-6, -3, 0, 3, 6), 
                  c("blue", "white", "white", "white", "red"))
col_freq <- colorRamp2(c(0, 1), c("white", "#e31a1c"))

## Set annotations
right_annotation <- rowAnnotation("Z-score" = anno_points(data_zscore,
                                                          ylim = c(0,5),
                                                          axis_param = list(side = "top",
                                                                            labels_rot = 0,
                                                                            gp = gpar(fontsize = 10))),
                                  "Coverage" = anno_points(data_coverage,
                                                           ylim = c(0,42),
                                                           axis_param = list(side = "top",
                                                                             labels_rot = 0,
                                                                             gp = gpar(fontsize = 10))),
                                  width = unit(3, "cm"),
                                  annotation_name_gp = gpar(fontsize = 13),
                                  annotation_name_side = "bottom",
                                  annotation_name_rot = 45,
                                  annotation_legend_param = list(border = TRUE))

## Set legend labels
heatmap_legend_param = list(title = "SD from healthy", 
                            border = TRUE,
                            at = c(-6, -3, 0, 3, 6))

## Set order
chr_order <- colnames(data_dist)
column_split <- data_chr$arm
row_order <- row.names(data_dist)
row_split <- data_cohort

## Generate oncopring
pdf(file.path(outdir, "Supplemental Figure 7.pdf"), height = 8, width = 14)
Heatmap <- Heatmap(data_dist,
                   col = col,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_row_dend = FALSE,
                   heatmap_legend_param = heatmap_legend_param,
                   right_annotation = right_annotation,
                   column_order = chr_order,
                   row_order = row_order,
                   row_labels = NULL,
                   column_split = column_split,
                   row_split = row_split, 
                   column_title_gp = gpar(fontsize = 10),
                   row_title_rot = 0,
                   border = TRUE,
                   border_gp = gpar(col = "black"))
draw(Heatmap, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

