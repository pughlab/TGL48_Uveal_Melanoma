library(tidyverse)
library(dplyr)
library(data.table)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)

### Set paths
path <- ""
outdir <- ""
PRC1 <- "H2AK119Ub_targets.txt"
healthy <- "HBC_fragment_ratios_100kb.txt"

frags <- file.path(path, "/fragmentomics/TGL48_UVM_fragment_ratio_100kb.txt")
samples <- file.path(path, "TGL48_UVM_sample_list.txt")

### Import data
data_frag <- read.delim(frags)
data_prc1 <- read.delim(PRC1)
samples <- read.delim(samples)
HBC <- read.delim(healthy, check.names = FALSE)

### Format sample list
samples <- samples[!(samples$Timepoint %in% c("Lymphocytes", "Tumour", "Healthy")), ]
samples$Patient <- gsub("UMB-0", "", samples$Patient)
samples <- samples[samples$Timepoint %in% c("Baseline", "2 weeks", "3 months", "6 months", "12 months"), ]
samples$Relapse <- factor(samples$Relapse, levels = c("No", "Yes", ""),
                          labels = c("No", "Yes", "LtFU"))
samples$timepoint <- as.character(samples$Timepoint)
samples$timepoint <- factor(samples$Timepoint,
                            levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))
samples <- samples[order(factor(samples$Relapse, levels = c("Yes", "No", "LOF")),
                         factor(samples$Stage, levels = c("IIB", "IIIA", "IIIB")),
                         factor(samples$Treatment, levels = c("Brachytherapy", "Enucleation")),
                         samples$Patient,
                         factor(samples$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))), ]

### Normalize and center 100kb bin ratios (Uveal cohort)
chr <- data_frag[ , c(1:3)]
data_frag <- as.matrix(data_frag[ , -c(1:3)])
data_frag[data_frag > 5] <- 5
data_frag[data_frag < -5] <- -5
data_frag <- ((data_frag - colMeans2(data_frag, na.rm = TRUE))/colSds(data_frag, na.rm = TRUE))*0.01
data_frag <- cbind(chr, data_frag)

### Normalize and center 100kb bin ratios (HBC cohort)
HBC <- as.matrix(HBC[ , -c(1:3)])
HBC <- scale(HBC)
HBC <- ((HBC - colMeans2(HBC, na.rm = TRUE))/colSds(HBC, na.rm = TRUE))*0.01
HBC <- cbind(chr, HBC)

### Find overlapping bins
colnames(data_prc1) <- c("ensemble", "gene", "seqnames", "start", "end")
prc1_targets <- as.data.table(data_prc1[ , c(3,4,5)])
frag_targets <- as.data.table(chr)

setDT(prc1_targets)
setkey(frag_targets)
overlap <- foverlaps(prc1_targets, frag_targets, type="within", nomatch=0L)

### Merge overlap with ratios
overlap <- overlap[ , c(1:3)]
overlap <- overlap[!duplicated(overlap),]
overlap$PRC1 <- "yes"
data_frag <- merge(data_frag, overlap, by = c("seqnames", "start", "end"), all = TRUE)
data_frag$PRC1[is.na(data_frag$PRC1)] <- "no"

HBC <- merge(HBC, overlap, by = c("seqnames", "start", "end"), all = TRUE)
HBC$PRC1[is.na(HBC$PRC1)] <- "no"

### Format ratio dataframe
data_frag <- data_frag[data_frag$PRC1 == "yes", ]
data_frag <- data_frag[ , samples$sWGS]
data_frag <- data_frag[complete.cases(data_frag), ]

### Make sure sample order is same as sample sheet
data_frag <- data_frag[ , colnames(data_frag) %in% samples$sWGS]
data_frag <- data_frag[ , samples$sWGS]

### Find healthy median of PRC1 Target Loci
HBC <- as.matrix(HBC[HBC$PRC1 == "yes", !(colnames(HBC) %in% c("seqnames", "start", "end", "PRC1"))])
HBC <- HBC[row.names(HBC) %in% row.names(data_frag), ]
HBC_median <- rowMedians(HBC, na.rm = TRUE)

### Calculate HBC Z-scores
mean <- mean(HBC_median)
sd <- sd(HBC_median)
HBC_means <- colMeans2(as.matrix(HBC), na.rm = TRUE)
HBC_sds <- colSds(as.matrix(HBC), na.rm = TRUE)
HBC_scores <- (HBC_means - mean)/sqrt((sd^2) + (HBC_sds)^2)
HBC_scores_mean <- mean(HBC_scores)
limit <- quantile(HBC_scores, 0.9)

### Calculate sample Z-scores
data_means <- colMeans2(as.matrix(data_frag), na.rm = TRUE)
data_sds <- colSds(as.matrix(data_frag), na.rm = TRUE)
data_scores <- (data_means - mean)/sqrt((sd^2) + (data_sds)^2)

### Attach HBC to Uveal cohort
data_frag$HBC <- HBC_median
data_frag <- as.matrix(data_frag)
HBC_sample <- c(replicate(5, "HBC"), NA, NA, replicate(5, "HBC"))
samples <- rbind(samples, HBC_sample)

### Set upper limits to cut off histogram
data_frag[data_frag > 0.1] <- 0.1
data_frag[data_frag < -0.1] <- -0.1

### Set Z-score annotation
data_scores <- as.data.frame(data_scores)
data_scores <- rbind(data_scores, HBC_scores_mean)
data_scores$limit <- limit
data_scores <- as.matrix(data_scores)
row.names(data_scores) <- colnames(data_frag)

### Set clinical information
data_relapse <- as.matrix(samples$Relapse)
row.names(data_relapse) <- samples$sWGS
data_relapse <- factor(data_relapse, levels = c("Yes", "No", ""),
                       labels = c("Relapse", "Remission", "Lost to Follow-up"))

data_time <- as.matrix(samples$Timepoint)
row.names(data_time) <- samples$TGL_ID
data_time <- factor(data_time, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))

### Set annotation and heatmap colours
col_fun <- colorRamp2(c(-0.05, 0, 0.05), 
                      c("#1f78b4", "white", "#e31a1c"))
col_relapse <- c(Relapse = "#fb9a99", Remission = "#a6cee3", "Lost to Follow-up" = "lightgrey", HBC = "lightgrey")
col_time <- c("Baseline" = "#FED976", "2 weeks" = "#FEB24C", "3 months" = "#FD8D3C", "6 months" = "#FC4E2A", "12 months" = "#E31A1C", HBC = "lightgrey")

### Create Annotation
PRC1_annotation <- HeatmapAnnotation("PRC1 Targets\nZ-Score" = anno_lines(data_scores, 
                                                            add_points = TRUE,
                                                            pch = c(16, NA),
                                                            gp = gpar (col = c("black", "red"),
                                                                       lty = c("solid", "dashed")),
                                                            pt_gp = gpar(col = c("black", NA)),
                                                            ylim = c(-1.25, 1.25),
                                                            axis_param = list(gp = gpar(fontsize = 10))),
                                    annotation_name_side = "left",
                                    annotation_name_gp = gpar(fontsize = 12),
                                    annotation_name_rot = 0,
                                    height = unit(1.75, "cm"),
                                    border = TRUE,
                                    annotation_legend_param = list(border = TRUE))

### Set column splits
col_split <- samples$Patient
split_order <- unique(col_split)
col_split <- factor(col_split, levels = split_order)
col_order <- colnames(data_frag)

### Set legend labels
heatmap_legend_param = list(title = "Fragment Ratio",
                            border = TRUE,
                            at = c(-0.05, 0, 0.05),
                            legend_height = unit(2, "cm"))

### Generate heatmap
pdf(file.path(outdir, "Figure 3 - Fragment ratio PRC1.pdf"), width = 12, height = 5)
PRC1_heatmap <- Heatmap(data_frag,
                        col = col_fun,
                        heatmap_legend_param = heatmap_legend_param,
                        column_order = col_order,
                        column_split = col_split,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        show_row_dend = FALSE,
                        top_annotation = PRC1_annotation,
                        border = FALSE,
                        heatmap_height = unit(3, "in"))
draw(PRC1_heatmap, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()

### Export z-scores
data_scores <- as.data.frame(data_scores)
data_scores <- data_scores[!(row.names(data_scores) == "HBC"), ]
data_scores$sample <- row.names(data_scores)
data_scores$limit <- (data_scores$limit - min(data_scores$data_scores))/(max(data_scores$data_scores) - min(data_scores$data_scores))
data_scores$PRC1 <- (data_scores$data_scores - min(data_scores$data_scores))/(max(data_scores$data_scores) - min(data_scores$data_scores))
data_scores <- data_scores[, c("sample", "PRC1", "limit")]
write.table(data_scores, file = file.path(outdir, "PRC1_scores.txt"), sep = "\t", row.names = FALSE)
