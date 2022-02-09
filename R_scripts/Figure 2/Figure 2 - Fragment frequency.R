library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project"
healthy_path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/combined_cohort"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Uveal Melanoma/Manuscript/Figures/Figure 2"

data_frequency <- file.path(path, "insert_size", "TGL48_UVM_fragment_swgs_freq.txt")
data_samples <- file.path(path, "TGL48_UVM_sample_list.txt")
data_healthy <- file.path(healthy_path, "insert_size", "HBC_fragment_freq.txt")

### Import data
data_frequency <- read.delim(data_frequency)
data_normal <- read.delim(data_healthy, check.names = FALSE)
data_samples <- read.delim(data_samples)

### Subset plasma samples and order based on clinical data 
data_samples <- data_samples[data_samples$sWGS %in% colnames(data_frequency), ]
data_samples$Patient <- gsub("UMB-0", "", data_samples$Patient)
data_samples$Timepoint <- factor(data_samples$Timepoint,
                                 levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))
data_samples <- data_samples[order(factor(data_samples$Relapse, levels = c("Yes", "No", ""),
                                          labels = c("Yes", "No", "Lost to Follow-up")),
                                   factor(data_samples$Stage, levels = c("IIB", "IIIA", "IIIB")),
                                   factor(data_samples$Treatment, levels = c("Brachytherapy", "Enucleation")),
                                   factor(data_samples$Volume),
                                   data_samples$Patient,
                                   factor(data_samples$Timepoint, levels = c("Tumour", "Baseline", "2 weeks", "3 months", "6 months", "12 months"))), ]

### Format fragment frequencies
row.names(data_frequency) <- data_frequency$length

### Make normal median
#outliers <- c("lib085", "lib150", "lib164")
#data_normal <- data_frequency[ , !(colnames(data_frequency) %in% data_samples$sWGS)]
#data_normal <- data_normal[ , !(colnames(data_normal) %in% outliers)]
row.names(data_normal) <- data_normal$length
data_normal <- data_normal[, -1]
data_normal <- as.matrix(data_normal)
data_frequency <- data_frequency[ , (colnames(data_frequency) %in% data_samples$sWGS)]
data_frequency <- data_frequency[ , data_samples$sWGS]
data_frequency <- as.matrix(data_frequency)
data_median <- rowMedians(data_normal)
data_sd <- rowSds(data_normal)
data_median <- as.data.frame(cbind(data_median, data_sd))
colnames(data_median) <- c("median", "sd")
row.names(data_median) <- row.names(data_frequency)
rm(data_sd)

### Normalize to fold change over median
data_change <- (data_frequency-data_median$median)/data_median$sd
data_change <- as.matrix(t(data_change))
data_change[is.nan(data_change)] <- 0
data_change[is.infinite(data_change)] <- 0

### Calculate summed Z-scores
data_normal <- (data_normal-data_median$median)/data_median$sd
healthy_sums <- colSums(abs(data_normal[c(1:141), ]), na.rm = TRUE)
healthy_median <- median(healthy_sums, na.rm = TRUE)
healthy_sd <- mad(healthy_sums)
data_sums <- rowSums(abs(data_change[ , c(1:141)]))
data_zscores <- as.data.frame((data_sums-healthy_median)/healthy_sd)

### Calculate Z-score threshold for 90% of healthy controls
healthy_zscores <- (healthy_sums - healthy_median)/healthy_sd
z_limit <- quantile(healthy_zscores, 0.9)

### Apphend limit and set Z-score annotation
data_zscores$limit <- z_limit
colnames(data_zscores) <- c("Z-score", "limit")
data_zscores <- as.matrix(data_zscores)

### Format heatmap table
data_change <- data_change[ , colnames(data_change) %in% c(10:320)]

### Set fragment annotations
data_fragment <- as.numeric(as.matrix(colnames(data_change)))
data_size <- as.matrix(ifelse(data_fragment <= 89, "1",
                  ifelse(data_fragment >=90 & data_fragment <=150, "2",
                         ifelse(data_fragment >=151 & data_fragment <=220, "3", "4"))))
row.names(data_size) <- colnames(data_change)
data_size <- factor(data_size,
                    levels = c("1", "2", "3", "4"),
                    labels = c("10-89", "90-150", "151-220", "221-320"))
rm(data_fragment)

### Set frequency distribution for the healthy median annotation
data_median <- data_median[row.names(data_median) %in% c(10:320), ]
data_median <- as.matrix(data_median$median)
row.names(data_median) <- colnames(data_change)

### Set clinical annotations
data_relapse <- as.matrix(data_samples$Relapse)
row.names(data_relapse) <- data_samples$sWGS
data_relapse <- factor(data_relapse, levels = c("Yes", "No", ""),
                       labels = c("Relapse", "Remission", "Lost to Follow-up"))

data_labels <- as.matrix(data_samples$Timepoint)
row.names(data_labels) <- data_samples$sWGS

### Set annotation and heatmap colours
col_fun <- colorRamp2(c(-6, -3, 0, 3, 6), 
                      c("#1f78b4", "white", "white", "white", "#e31a1c"))
col_relapse <- c(Relapse = "#fb9a99", Remission = "#a6cee3", "Lost to Follow-up" = "lightgrey")

### Compile annotation layers
left_annotation <- rowAnnotation(Outcome = data_relapse,
                                 show_annotation_name = TRUE,
                                 annotation_name_side = "top",
                                 border = TRUE,
                                 col = list(Outcome = col_relapse),
                                 annotation_legend_param = list(border = TRUE,
                                                                legend_direction = "horizontal",
                                                                nrow = 1),
                                 simple_anno_size = unit(0.3, "cm"))

right_annotation <- rowAnnotation("Z-Score" = anno_lines(data_zscores,
                                                         add_points = TRUE,
                                                         pch = c(16, NA),
                                                         gp = gpar(col = c("white", "red"),
                                                         lty = c("solid", "dashed")),
                                                         ylim = c(0, 17),
                                                         axis_param = list(side = "top",
                                                                           labels_rot = 0,
                                                                           gp = gpar(fontsize = 10)),
                                                         width = unit(1.5, "cm")),
                                  Timepoint = anno_text(data_labels,
                                                        gp = gpar(fontsize = 8)),
                                  border = TRUE,
                                  annotation_name_side = "top",
                                  annotation_name_rot = 0,
                                  width = unit(4, "cm"))

freq_annotation <- HeatmapAnnotation("Healthy Median" = anno_lines(data_median,
                                                                   gp = gpar(fontsize = 10),
                                                                   axis = FALSE),
                                     border = TRUE,
                                     annotation_name_side = "right",
                                     annotation_name_gp = gpar(fontsize = 14),
                                     height = unit(3, "cm"),
                                     annotation_name_rot = 0)

### Set legend labels
heatmap_legend_param = list(title = "SD from Healthy",
                            border = TRUE,
                            at = c(-6, -3, 0, 3, 6),
                            direction = "horizontal",
                            legend_width = unit(3, "cm"))

### Set column/row orders and splits
col_order <- colnames(data_change)
row_order <- data_samples$sWGS
row_split <- data_samples$Patient
split_order <- unique(row_split)
row_split <- factor(row_split, levels = split_order)
col_split <- data_size

## Generate heatmap
pdf(file.path(outdir, "Figure 2 - Fragment frequency.pdf"), width = 5, height = 7)
Heatmap <- Heatmap(data_change,
                   col = col_fun,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_order = row_order,
                   column_order = col_order,
                   left_annotation = left_annotation,
                   right_annotation = right_annotation,
                   top_annotation = freq_annotation,
                   heatmap_legend_param = heatmap_legend_param,
                   row_split = row_split,
                   column_split = col_split,
                   row_title_rot = 0,
                   row_title_gp = gpar(fontsize = 10),
                   column_title_gp = gpar(fontsize = 10),
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   border = TRUE)
draw(Heatmap, column_title = "Fragment Size (bp)", row_title = "Patient", heatmap_legend_side = "bottom", merge_legend = TRUE) 
dev.off()

### Prepare metrics for export
data_z <- as.data.frame(data_zscores)
data_z$sample <- row.names(data_z)
data_z$limit <- (data_z$limit - min(data_z$`Z-score`))/(max(data_z$`Z-score`) - min(data_z$`Z-score`))
data_z$insert <- (data_z$`Z-score` - min(data_z$`Z-score`))/(max(data_z$`Z-score`) - min(data_z$`Z-score`))
data_z <- data_z[, c("sample", "insert", "limit")]
write.table(data_z, file.path(outdir, "insert_size.txt"), sep = "\t", row.names = FALSE)
