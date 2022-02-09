library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(matrixStats)

### Set paths
sites <- c("Liver_dnase")
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Uveal Melanoma/Manuscript/Figures/Figure 4"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/combined_cohort/griffin/uveal"

data_griffin <- file.path(path, "griffin/uveal", paste0("TGL48_UVM_griffin_corrected_", sites, ".txt"))
data_samples <- file.path(path, "TGL48_UVM_sample_list.txt")
data_normal <- file.path(healthy_path, paste0("HBC_griffin_corrected_", sites, ".txt"))

### Import data 
data_griffin <- lapply(data_griffin, read.delim)
data_normal <- lapply(data_normal, read.delim)
data_samples <- read.delim(data_samples)

### Keep only plasma samples and order based on clinical information
data_samples <- data_samples[!(data_samples$Timepoint %in% c("Lymphocytes", "Tumour", "Healthy")), ]
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

### Name dataframes in data list
names(data_griffin) <- sites
names(data_normal) <- paste0("normal_", sites)

### Format and order dataframes
my_fun <- function(x){
row.names(x) <- x$distance
x <- x[, -c(1)]
x <- x[, data_samples$sWGS]
x <- as.matrix(x)
}
data_griffin <- lapply(data_griffin, my_fun)

my_fun <- function(x){
row.names(x) <- x$distance
x <- as.matrix(x[, -c(1)])
}
data_normal <- lapply(data_normal, my_fun)

### Center data
my_fun <- function(x) {
  x <- scale(x)*0.1 + 1
}
data_griffin <- lapply(data_griffin, my_fun)
data_normal <- lapply(data_normal, my_fun)

### Calculate the healthy median
normal_median <- lapply(data_normal, rowMedians)
names(normal_median) <- paste0(sites, "_healthy_median")

normal_sd <- lapply(data_normal, rowMads)
names(normal_sd) <- paste0(sites, "_healthy_sd")

### Split dataframes into each site
list2env(data_griffin, envir = .GlobalEnv)
list2env(data_normal, envir = .GlobalEnv)
list2env(normal_median, envir = .GlobalEnv)
list2env(normal_sd, envir = .GlobalEnv)

### Calculate the normal mean (+/- 30bp) and healthy z-scores
normal_Liver_dnase_means <- colMeans2(normal_Liver_dnase[row.names(Liver_dnase) %in% c("-30", "-15", "0", "15", "30"), ])
normal_Liver_dnase_mean <- mean(normal_Liver_dnase_means)
normal_Liver_dnase_sd <- sd(normal_Liver_dnase_means)
normal_Liver_dnase_scores <- -(normal_Liver_dnase_means - normal_Liver_dnase_mean)/normal_Liver_dnase_sd
normal_Liver_dnase_limit <- quantile(normal_Liver_dnase_scores, 0.9)

### Calculate sample means (+/- 30bp) and z-scores
Liver_dnase_means <- colMeans2(Liver_dnase[row.names(Liver_dnase) %in% c("-30", "-15", "0", "15", "30"), ])
Liver_dnase_zscore <- -(Liver_dnase_means - normal_Liver_dnase_mean)/normal_Liver_dnase_sd

### Calculate the distance from the healthy median
Liver_dnase <- -(Liver_dnase - Liver_dnase_healthy_median)/Liver_dnase_healthy_sd

### Set healthy medians 
Liver_dnase_healthy <- as.matrix(Liver_dnase_healthy_median)
row.names(Liver_dnase_healthy) <- row.names(Liver_dnase)

### Set clinical information
data_relapse <- as.matrix(data_samples$Relapse)
row.names(data_relapse) <- row.names(data_samples)
data_relapse <- factor(data_relapse, levels = c("Yes", "No", ""),
                       labels = c("Relapse", "Remission", "Lost to Follow-up"))

data_time <- as.matrix(data_samples$Timepoint)
row.names(data_time) <- data_samples$TGL_ID
data_time <- factor(data_time, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))

### Add unique rownames to each sample table
row.names(Liver_dnase) <- paste0("dnase", rownames(Liver_dnase))

### Concatenate samples and format for Heatmap
data_griffin <- Liver_dnase
data_griffin <- t(data_griffin)

cohort_sites <- as.matrix(c(rep("Liver", nrow(Liver_dnase))))
row.names(cohort_sites) <- colnames(data_griffin)

### Concatenate healthy medians and format for Heatmap
data_healthy <- Liver_dnase_healthy

### Concatenate zscores
liver_limit <- c(rep(normal_Liver_dnase_limit, nrow(data_samples)))
data_zscore <- as.matrix(cbind(Liver_dnase_zscore, liver_limit))

### Focus on +/- 30bp of the site
data_griffin[, !(colnames(data_griffin) %in% c("dnase-30", "dnase-15", "dnase0", "dnase15", "dnase30"))] <- 0

### Set annotation and heatmap colours
col <- colorRamp2(c(-3, 0, 3), 
                  c("#1f78b4", "white", "#e31a1c"))
col_relapse <- c(Relapse = "#fb9a99", Remission = "#a6cee3", "Lost to Follow-up" = "lightgrey")
col_time <- c("Baseline" = "#FED976", "2 weeks" = "#FEB24C", "3 months" = "#FD8D3C", "6 months" = "#FC4E2A", "12 months" = "#E31A1C")

### Set additional annotations
left_annotation <- rowAnnotation("Outcome" = data_relapse,
                                 "Timepoint" = data_time,
                                 show_annotation_name = FALSE,
                                 border = TRUE,
                                 col = list("Outcome" = col_relapse, "Timepoint" = col_time),
                                 annotation_legend_param = list(border = TRUE,
                                                                legend_direction = "vertical",
                                                                legend_label_gp = gpar(fontsize = 1)),
                                 simple_anno_size = unit(0.3, "cm"))

right_annotation <- rowAnnotation("Score" = anno_lines(data_zscore,
                                                       add_points = TRUE,
                                                       pch = c(16, NA),
                                                       gp = gpar(col = c("black", "red"),
                                                                lty = c("solid", "dashed")),
                                                       pt_gp = gpar(col = c("black", "red"),
                                                                    size = c(2, 0)),
                                                       ylim = c(-2, 3),
                                                       size = unit(1, "mm"),
                                                       axis_param = list(side = "top",
                                                                         labels_rot = 0,
                                                                         gp = gpar(fontsize = 8))),
                                  width = unit(2, "cm"),
                                  border = FALSE,
                                  annotation_name_gp = gpar(fontsize = 10),
                                  annotation_name_side = "top",
                                  annotation_name_rot = 0,
                                  annotation_legend_param = list(border = TRUE,
                                                                 legend_direction = "vertical"))

top_annotation <- HeatmapAnnotation("Healthy\nMedian" = anno_lines(data_healthy,
                                                                   ylim = c(0.8, 1.2),
                                                                   axis_param = list(side = "left",
                                                                                     labels_rot = 0,
                                                                                     gp = gpar(fontsize = 8)),
                                                                  border = FALSE),
                                    height = unit(1, "cm"),
                                    annotation_name_side = "left",
                                    annotation_name_rot = 0,
                                    annotation_name_gp = gpar(fontsize = 10))

## Set legend labels
heatmap_legend_param = list(title = "SD from healthy", 
                            legend_gp = gpar(border = "black"),
                            direction = "vertical",
                            border = TRUE)

## Set column/row order and column/row splits
sample_order <- data_samples$sWGS
bin_order <- colnames(data_griffin)

row_split <- data_samples$Patient
split_order <- unique(row_split)
row_split <- factor(row_split, levels = split_order)

col_split <- cohort_sites
split_order <- unique(col_split)
col_split <- factor(col_split, levels = split_order)

## Generate oncoprint
pdf(file.path(outdir, "Figure 4 - Griffin uveal sites.pdf"), height = 8, width = 4)
Griffin <- Heatmap(data_griffin,
                   col = col,
                   show_row_names = FALSE, 
                   show_column_names = FALSE,
                   heatmap_legend_param = heatmap_legend_param,
                   left_annotation = left_annotation,
                   right_annotation = right_annotation,
                   top_annotation = top_annotation,
                   row_order = sample_order,
                   column_order = bin_order,
                   row_labels = NULL,
                   row_title_gp = gpar(fontsize = 10),
                   column_title_gp = gpar(fontsize = 10),
                   row_split = row_split,
                   column_split = cohort_sites,
                   row_title_rot = 0,
                   border = TRUE)
draw(Griffin, heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()

### Prepare metrics for export
data_z <- as.data.frame(data_zscore)
data_z$sample <- data_samples$sWGS

data_z$liver_limit <- (data_z$liver_limit - min(data_z$Liver_dnase_zscore))/(max(data_z$Liver_dnase_zscore) - min(data_z$Liver_dnase_zscore))
data_z$liver <- (data_z$Liver_dnase_zscore - min(data_z$Liver_dnase_zscore))/(max(data_z$Liver_dnase_zscore) - min(data_z$Liver_dnase_zscore))
data_z <- data_z[, c("sample",  "liver", "liver_limit")]
write.table(data_z, file.path(outdir, "griffin_metrics.txt"), sep = "\t")
