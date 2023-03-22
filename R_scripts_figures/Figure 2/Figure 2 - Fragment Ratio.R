library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(matrixStats)

### Set variables
#path <- ""
#healthy_path <- ""
#outdir <- ""

data_frequency <- file.path("Input_data", "TGL48_UVM_fragment_ratio_5Mb.txt")
data_samples <- file.path("Input_data", "TGL48_UVM_sample_list.txt")
data_healthy <- file.path("Healthy_controls", "HBC_fragment_ratios.txt")

### Import data (Starting with the 5Mb ratios)
data_ratio <- read.delim(data_frequency)
data_normal <- read.delim(data_healthy, check.names = FALSE)
data_samples <- read.delim(data_samples)

### Keep only plasma samples and order based on clinical information
data_samples$Patient <- gsub("UMB-0", "", data_samples$Patient)
data_samples <- data_samples[data_samples$Timepoint %in% c("Baseline", "2 weeks", "3 months", "6 months", "12 months"), ]
data_samples$timepoint <- as.character(data_samples$Timepoint)
data_samples$timepoint <- factor(data_samples$Timepoint,
                                 levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))
data_samples <- data_samples[order(factor(data_samples$Relapse, levels = c("Yes", "No", "LOF")),
                                   factor(data_samples$Stage, levels = c("IIB", "IIIA", "IIIB")),
                                   factor(data_samples$Treatment, levels = c("Brachytherapy", "Enucleation")),
                                   data_samples$Patient,
                                   factor(data_samples$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))), ]

### Set chromosomes
data_chr <- data_ratio[, c("seqnames", "arm", "start", "end")]
row.names(data_chr) <- with(data_chr, paste0(seqnames, "_", start))

### Format ratio sheet and order samples properly
row.names(data_ratio) <- with(data_ratio, paste0(seqnames, "_", start))
data_ratio <- data_ratio[, -c(1:4)]
data_ratio <- as.matrix(data_ratio[ , colnames(data_ratio) %in% data_samples$sWGS])
data_ratio <- data_ratio[ , data_samples$sWGS]

### Calculate the healthy median
row.names(data_normal) <- with(data_normal, paste0(seqnames, "_", start))
data_normal <- as.matrix(data_normal[, -c(1:4)])
normal_median <- rowMedians(data_normal, na.rm = TRUE)
normal_sd <- rowSds(data_normal, na.rm = TRUE)

### Calculate the distance from the healthy median
data_ratio <- (data_ratio - normal_median)/normal_sd
data_normal <- (data_normal - normal_median)/normal_sd

### Calculate normal Z-scores
normal_sum <- colSums(abs(data_normal), na.rm = TRUE)
normal_sum_median <- median(normal_sum, na.rm = TRUE)
normal_MAD <- mad(normal_sum, na.rm = TRUE)

### Calculate Z-scores and set annotation
healthy_zscores <- abs((normal_sum - normal_sum_median))/normal_MAD
healthy_mean <- mean(healthy_zscores)
z_limit <- quantile(healthy_zscores, 0.9)

data_sum <- colSums(abs(data_ratio), na.rm = TRUE)
data_zscore <- (data_sum - normal_sum_median)/normal_MAD

data_zscore <- as.data.frame(data_zscore)
data_zscore <- rbind(data_zscore, healthy_mean) ### This is for the healthy control
data_zscore$limit <- z_limit
data_zscore <- as.matrix(data_zscore)

### Set healthy median 
lower <- normal_median - normal_sd
upper <- normal_median + normal_sd
normal_median <- as.matrix(cbind(normal_median, lower, upper))
row.names(normal_median) <- row.names(data_ratio)

### Set alteration frequencies
data_freq <- rowMeans2(data_ratio, na.rm = TRUE)

### Attach HBC to Uveal cohort
data_ratio <- as.data.frame(data_ratio)
data_ratio$HBC <- 0
data_ratio <- as.matrix(data_ratio)
HBC_sample <- c(replicate(5, "HBC"), NA, NA, NA, replicate(3, "HBC"), replicate(6, NA), "HBC")
data_samples <- rbind(data_samples, HBC_sample)

### Set clinical information
data_relapse <- as.matrix(data_samples$Relapse)
row.names(data_relapse) <- row.names(data_samples)
data_relapse <- factor(data_relapse, levels = c("Yes", "No", ""),
                       labels = c("Relapse", "Remission", "Lost to Follow-up"))

data_time <- as.matrix(data_samples$Timepoint)
row.names(data_time) <- data_samples$TGL_ID
data_time <- factor(data_time, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))

### Set annotation and heatmap colours
col_heat <- colorRamp2(c(-6, -1.5, 0, 1.5, 6), 
                      c("#1f78b4", "white", "white", "white", "#e31a1c"))
col_fun <- colorRamp2(c(-2, 0, 2), 
                       c("#1f78b4", "white", "#e31a1c"))
col_relapse <- c(Relapse = "#fb9a99", Remission = "#a6cee3", "Lost to Follow-up" = "lightgrey")
col_time <- c("Baseline" = "#FED976", "2 weeks" = "#FEB24C", "3 months" = "#FD8D3C", "6 months" = "#FC4E2A", "12 months" = "#E31A1C")

### Set additional annotations
top_annotation <- HeatmapAnnotation("Outcome" = data_relapse,
                                    "Timepoint" = data_time,
                                    show_annotation_name = TRUE,
                                    border = TRUE,
                                    col = list("Outcome" = col_relapse, "Timepoint" = col_time),
                                    annotation_name_gp = gpar(fontsize = 12),
                                    annotation_name_side = "right",
                                    annotation_name_rot = 0,
                                    annotation_legend_param = list(border = TRUE),
                                    simple_anno_size = unit(0.35, "cm"))

bot_annotation <- HeatmapAnnotation("Genome-Wide\nZ-score" = anno_lines(data_zscore,
                                                                        add_points = TRUE,
                                                                        pch = c(16, NA),
                                                                        gp = gpar(col = c("black", "red"),
                                                                                  lty = c("solid", "dashed")),
                                                                        ylim = c(0, 7),
                                                                        axis_param = list(side = "left",
                                                                                          labels_rot = 0,
                                                                                          gp = gpar(fontsize = 10))),
                                    height = unit(1.75, "cm"),
                                    show_annotation_name = TRUE,
                                    border = TRUE,
                                    annotation_name_gp = gpar(fontsize = 12),
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    annotation_legend_param = list(border = TRUE))

freq_annotation <- rowAnnotation("Cohort Mean" = data_freq,
                                 " " = anno_mark(at = c(84, 104, 207, 286, 288, 293, 299, 452, 480),
                                                 labels = c("SF3B1", "BAP1", "MAPK14", "MYC", "PTK2", "CDKN2A/B", "GNAQ", "TP53", "GNA11"),
                                                 which = "row",
                                                 side = "right",
                                                 labels_rot = 0,
                                                 padding = unit(5, "mm"),
                                                 labels_gp = gpar(fontsize = 10)),
                                 col = list("Cohort Mean" = col_fun),
                                 border = FALSE,
                                 annotation_name_side = "bottom",
                                 annotation_name_rot = 90,
                                 annotation_name_gp = gpar(fontsize = 13),
                                 annotation_legend_param = list(border = TRUE,
                                                                legend_height = unit(2.5, "cm")))

left_annotation <- rowAnnotation("Healthy\nMedian" = anno_lines(normal_median,
                                                                ylim = c(-0.08, 0.08),
                                                                axis_param = list(at = c(-0.05, 0, 0.05),
                                                                                  side = "top",
                                                                                  labels_rot = 0,
                                                                                  gp = gpar(fontsize = 8)),
                                                                border = FALSE),
                                 width = unit(2, "cm"),
                                 show_annotation_name = FALSE,
                                 annotation_name_rot = 0,
                                 annotation_name_gp = gpar(fontsize = 13))

### Set legend labels
heatmap_legend_param <- list(title = "SD from Healthy",
                             at = c(-6, -3, 0, 3, 6),
                             legend_gp = gpar(border = "black"),
                             border = TRUE,
                             legend_height = unit(2.5, "cm"))

## Set chromosome order
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

## Set column/row order and column/row splits
sample_order <- row.names(data_samples)
chr_order <- rownames(data_ratio)
column_split <- data_samples$Patient
split_order <- unique(column_split)
column_split <- factor(column_split, levels = split_order)
row_split <- data_chr$arm

## Generate oncoprint
pdf(file.path(outdir, "Figure 3 - Fragment ratio.pdf"), height = 8, width = 10)
Heatmap <- Heatmap(data_ratio,
                   col = col_heat,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   heatmap_legend_param = heatmap_legend_param,
                   top_annotation = top_annotation,
                   bottom_annotation = bot_annotation,
                   left_annotation = left_annotation,
                   right_annotation = freq_annotation,
                   row_order = chr_order,
                   column_order = data_samples$sWGS,
                   row_labels = NULL,
                   column_split = column_split,
                   column_title_gp = gpar(fontsize = 10),
                   row_title_gp = gpar(fontsize = 8),
                   row_split = row_split,
                   row_title_rot = 0,
                   border = FALSE,
                   heatmap_height = unit(6, "in"))
draw(Heatmap, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

### Prepare metrics for export
data_z <- as.data.frame(data_zscore)
data_z <- data_z[rownames(data_z) %in% data_samples$sWGS, ]
data_z$sample <- row.names(data_z)
data_z$limit <- (data_z$limit - min(data_z$data_zscore))/(max(data_z$data_zscore) - min(data_z$data_zscore))
data_z$ratio <- (data_z$data_zscore - min(data_z$data_zscore))/(max(data_z$data_zscore) - min(data_z$data_zscore))
data_z <- data_z[, c("sample", "ratio", "limit")]
write.table(data_z, file.path(outdir, "ratios.txt"), sep = "\t", row.names = FALSE)
