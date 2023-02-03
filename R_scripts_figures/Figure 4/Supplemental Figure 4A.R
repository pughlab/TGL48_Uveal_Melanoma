rm(list = ls())
cat("\014")

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(tidyverse))
library(dplyr)
library(data.table)
suppressPackageStartupMessages(library(circlize))

setwd('')

# medip data
df <- read_tsv('inputdata/uveal.tsv', show_col_types = F)
df[,c("bin_chr","bin_start","bin_end","i.bin_start","i.bin_end")] <- NULL
data <- df[,47:53]
df <- read.delim('inputdata/metadata.tsv', colClasses=c("character", "factor", "factor"))
metadata <- df[47:53,]
score <- colSums(data)
score <- as.matrix(score)
data <- data[, order(factor(metadata$Patient,levels=c("09","03","05","04","06","01","02","11","07","08","10")))]
score <- score[order(factor(metadata$Patient,levels=c("09","03","05","04","06","01","02","11","07","08","10"))),]
metadata <- metadata[order(factor(metadata$Patient,levels=c("09","03","05","04","06","01","02","11","07","08","10"))),]
data <- as.matrix(data)
score <- as.matrix(score)

# new matrix --- uveal, HBC, Kui
df_new <- read_tsv('inputdata/uveal.tsv', show_col_types = F)
df_new[,c("bin_chr","bin_start","bin_end","i.bin_start","i.bin_end")] <- NULL
HBC <- df_new[,54:67]
# data <- full_join(data, HBC, by=c('bin_chr', 'bin_start', 'bin_end'))

# patientHBC <- c('HBC-01-0015','HBC-01-0011','HBC-01-0014','HBC-01-0013','HBC-01-0012','HBC-01-0010','HBC-01-0009',
#                 'HBC-01-0017','HBC-01-0001','HBC-01-0002','HBC-01-0003','HBC-01-0004','HBC-01-0005','HBC-01-0006')
patientHBC <- c('15','11','14','13','12','10','09','17','01','02','03','04','05','06')
timepointHBC <- rep('HBC', 14)
relapseHBC <- rep('HBC', 14)
metaHBC <- data.frame(patientHBC,timepointHBC,relapseHBC)
colnames(metaHBC) <- colnames(metadata)

scoreHBC <- colSums(HBC)
scoreHBC <- as.matrix(scoreHBC)
HBC <- HBC[, order(factor(metaHBC$Patient,levels=c('01','02','03','04','05','06','09','10','11','12','13','14','15','17')))]
scoreHBC <- scoreHBC[order(factor(metaHBC$Patient,levels=c('01','02','03','04','05','06','09','10','11','12','13','14','15','17'))),]
metaHBC <- metaHBC[order(factor(metaHBC$Patient,levels=c('01','02','03','04','05','06','09','10','11','12','13','14','15','17'))),]
HBC <- as.matrix(HBC)
scoreHBC <- as.matrix(scoreHBC)

# combine together
# data <- cbind(data, HBC)
# score <- rbind(score, scoreHBC)
# metadata <- rbind(metadata, metaHBC)

# score to derek
# HBC.median <- median(scoreHBC)
# HBC.sd <- sqrt(var(scoreHBC))
# tumor.median <- median(score)
# tumor.sd <- sqrt(var(score))
# write.table(score,'med2/tumor_score_0527.tsv',quote=F, sep='\t')
################

############--------------------------------------############
data_relapse <- as.matrix(metadata$Relapse)
row.names(data_relapse) <- colnames(data)
data_relapse <- factor(data_relapse, levels = c("Yes", "No", "LOF", "HBC"),
                       labels = c("Relapse", "Remission", "Lost to Follow-up", "HBC"))
col_relapse <- c(Relapse = "#fb9a99", Remission = "#a6cee3", "Lost to Follow-up" = "lightgrey",
                 HBC = "#b2df8a")

top_ha <- HeatmapAnnotation(Outcome = data_relapse,
                            "Methylation\nScore" = anno_lines(score, add_points = TRUE,
                                                              ylim=c(10,70),
                                                              pch = c(16),
                                                              gp = gpar (col = c("black"),
                                                                         lty = c("solid")),
                                                              pt_gp = gpar(col = c("black")),
                                                              axis_param = list(gp = gpar(fontsize = 10))),
                            col = list(Outcome = col_relapse),
                            show_annotation_name = c("Methylation\nScore" = FALSE, Outcome = FALSE),
                            annotation_name_side = "right",
                            annotation_name_gp = gpar(fontsize = 12),
                            annotation_name_rot = 0,
                            height = unit(3, "cm"),
                            border = TRUE,
                            annotation_legend_param = list(border = TRUE))

col_split <- metadata$Patient
split_order <- unique(col_split)
col_split <- factor(col_split, levels = split_order)
col_order <- colnames(data)

heatmap_legend_param = list(title = "Methylation Level",
                            border = TRUE,
                            at = c(0, 0.5, 1),
                            legend_height = unit(2, "cm"))

col_fun <- colorRamp2(c(0, 1), c("white", "#e31a1c"))
pdf("suppl_Figure_4A.pdf", width = 7, height = 7)
cfmedip_heatmap <- Heatmap(data,
                           col = col_fun,
                           heatmap_legend_param = heatmap_legend_param,
                           column_order = as.vector(metadata$V1),
                           column_split = col_split,
                           show_row_names = FALSE,
                           show_column_names = FALSE,
                           show_row_dend = FALSE,
                           top_annotation = top_ha,
                           border = FALSE,
                           heatmap_height = unit(5, "in"))




# HBC heatmap
data_relapseHBC <- as.matrix(metaHBC$Relapse)
row.names(data_relapseHBC) <- colnames(HBC)
data_relapseHBC <- factor(data_relapseHBC, levels = c("Yes", "No", "LOF", "HBC"),
                       labels = c("Relapse", "Remission", "Lost to Follow-up", "HBC"))
col_relapseHBC <- c(Relapse = "#fb9a99", Remission = "#a6cee3", "Lost to Follow-up" = "lightgrey",
                 HBC = "#b2df8a")

top_haHBC <- HeatmapAnnotation(Outcome = data_relapseHBC,
                            "Methylation\nScore" = anno_lines(scoreHBC, add_points = TRUE,
                                                              ylim=c(10,70),
                                                              pch = c(16),
                                                              gp = gpar (col = c("black"),
                                                                         lty = c("solid")),
                                                              pt_gp = gpar(col = c("black")),
                                                              axis_param = list(gp = gpar(fontsize = 10), side='right')),
                            col = list(Outcome = col_relapseHBC),
                            annotation_name_side = "right",
                            annotation_name_gp = gpar(fontsize = 12),
                            annotation_name_rot = 0,
                            height = unit(3, "cm"),
                            border = TRUE,
                            show_legend = F,
                            annotation_legend_param = list(border = TRUE))

col_splitHBC <- metaHBC$Patient
split_orderHBC <- unique(col_splitHBC)
col_splitHBC <- factor(col_splitHBC, levels = split_orderHBC)
col_orderHBC <- colnames(HBC)

heatmap_legend_param = list(title = "Methylation Level",
                            border = TRUE,
                            at = c(0, 0.5, 1),
                            legend_height = unit(2, "cm"))

col_fun <- colorRamp2(c(0, 1), c("white", "#e31a1c"))

HBC_heatmap <- Heatmap(HBC,
                       col = col_fun,
                       show_heatmap_legend = FALSE,
                       column_order = as.vector(metaHBC$V1),
                           column_split = col_splitHBC,
                           show_row_names = FALSE,
                           show_column_names = FALSE,
                           show_row_dend = FALSE,
                           top_annotation = top_haHBC,
                           border = FALSE,
                           heatmap_height = unit(5, "in"))
result <- cfmedip_heatmap + HBC_heatmap
draw(result, merge_legend = TRUE, heatmap_legend_side = "right", ht_gap = unit(0.5, "cm"))
dev.off()
