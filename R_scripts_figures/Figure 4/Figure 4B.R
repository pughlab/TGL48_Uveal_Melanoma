rm(list = ls())
cat("\014")

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(tidyverse))
library(dplyr)
library(data.table)
suppressPackageStartupMessages(library(circlize))

setwd('/Users/pluo/Project/UVM/')

data <- read_tsv('inputdata/uveal_liver_sort.tsv', show_col_types = F)

df <- read.delim('inputdata/metadata.tsv', colClasses=c("character", "factor", "factor"))
metadata <- df[1:46,]

# HBC metadata
temp <- data.frame('HBC','HBC','HBC')
colnames(temp) <- colnames(metadata)
metadata <- rbind(metadata, temp)

# HBC score
df_new <- read_tsv('inputdata/uveal_liver.tsv', show_col_types = F)
df_new[,c("bin_chr","bin_start","bin_end","i.bin_start","i.bin_end")] <- NULL
HBC <- df_new[,54:67]
scoreHBC <- colSums(HBC)
medianHBC <- median(scoreHBC)
medianHBC <- as.data.frame(medianHBC)
rownames(medianHBC) <- "HBC"
colnames(medianHBC) <- "score"
HBC90 <- quantile(scoreHBC, c(.90))

score <- colSums(data)
score <- as.data.frame(score)
score <- rbind(score, medianHBC)

#--------# to Derek
# Derek_score <- cbind(score, metadata)
# write.table(Derek_score, 'med2/liver_score_0527.tsv', quote=F, sep='\t')
#--------#

score$baseline <- rep(HBC90, 47)
# score$baseline <- c(rep(score[1,],4),rep(score[5,],4),rep(score[9,],4),rep(score[13,],5),rep(score[18,],4),
#                     rep(score[22,],5),rep(score[27,],5),rep(score[32,],5),rep(score[37,],4),score[41,],
#                     rep(score[42,],5))
score <- as.matrix(score)

HBC <- rowMeans(HBC)
data <- cbind(data, HBC)

# sort data matrix
data <- data[, order(factor(metadata$Patient,levels=c("09","03","05","04","06","01","02","11","07","08","10","HBC")))]
score <- score[order(factor(metadata$Patient,levels=c("09","03","05","04","06","01","02","11","07","08","10","HBC"))),]
metadata <- metadata[order(factor(metadata$Patient,levels=c("09","03","05","04","06","01","02","11","07","08","10","HBC"))),]

data <- as.matrix(data)
score <- as.matrix(score)

data_relapse <- as.matrix(metadata$Relapse)
row.names(data_relapse) <- colnames(data)
data_relapse <- factor(data_relapse, levels = c("Yes", "No", "LOF", "HBC"),
                       labels = c("Relapse", "Remission", "Lost to Follow-up", "HBC"))

data_time <- as.matrix(metadata$Timepoint)
row.names(data_time) <- colnames(data)
data_time <- factor(data_time, levels = c("Baseline","2 weeks","3 months","6 months","12 months", "HBC"))

col_relapse <- c(Relapse = "#fb9a99", Remission = "#a6cee3", "Lost to Follow-up" = "lightgrey", "HBC" = "lightgrey")
col_time <- c("Baseline" = "#FED976", "2 weeks" = "#FEB24C", "3 months" = "#FD8D3C", "6 months" = "#FC4E2A",
              "12 months" = "#E31A1C", "HBC" = "lightgrey")

top_ha <- HeatmapAnnotation(Outcome = data_relapse,
                            Timepoint = data_time,
                            "Methylation\nScore" = anno_lines(score, add_points = TRUE,
                                                              pch = c(16, NA),
                                                              gp = gpar (col = c("black","red"),
                                                                     lty = c("solid","dashed")),
                                                              pt_gp = gpar(col = c("black",NA)),
#                                                             ylim = c(-1.25, 1.25),
                                                          axis_param = list(gp = gpar(fontsize = 10))),
                            col = list(Outcome = col_relapse, Timepoint = col_time),
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
pdf("Figure_4B.pdf", width = 12, height = 7)
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
draw(cfmedip_heatmap, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()

#HBC.median <- median(scoreHBC)
#HBC.sd <- sqrt(var(scoreHBC))
#tumor.median <- median(score)
#tumor.sd <- sqrt(var(score))

#remission median and sd
# score_remission <- score[metadata$Relapse == 'No',]
# score_remission <- as.data.frame(score_remission)
# median(score_remission$score)
# sqrt(var(score_remission$score))
