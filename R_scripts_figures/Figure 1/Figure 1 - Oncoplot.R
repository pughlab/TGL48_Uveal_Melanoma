library(ComplexHeatmap)
library(dplyr)
library(matrixStats)
library(circlize)

### Set variables
path <- ""
healthy_path <- ""
outdir <- ""

data_onco <- file.path(path, "pipeline_analysis", "combined_variants.txt")
data_samples <- file.path(path, "TGL48_UVM_sample_list.txt")
data_depths <- file.path(path, "ichorCNA", "TGL48_UVM_corrected_depths.txt")
data_depths_short <- file.path(path, "ichorCNA", "TGL48_UVM_corrected_depths_short.txt")
healthy_depths <- file.path(healthy_path, "ichorCNA", "HBC_ichorCNA_corrected_depths.txt")
healthy_depths_short <- file.path(healthy_path, "ichorCNA", "HBC_ichorCNA_corrected_depths_short.txt")

### Import data
data_onco <- read.delim(data_onco)
data_samples <- read.delim(data_samples)
data_depths <- read.delim(data_depths)
data_depths_short <- read.delim(data_depths_short)
healthy_depths <- read.delim(healthy_depths)
healthy_depths_short <- read.delim(healthy_depths_short)

### Generate Z-score statistics for ichorCNA data
### Remove problem bins and NA bins
data_depths <- merge(data_depths, healthy_depths, by = c("chr", "start", "end"), all = TRUE)
problems <- c(285, 434, 677, 710, 711, 983, 984, 1187, 1263) ## These bins have abnormal drops and skew data
data_depths <- data_depths[!(rownames(data_depths) %in% problems), -c(2,3)]
data_depths <- data_depths[complete.cases(data_depths), ]

### Normalize copy number across chromosomes
names <- colnames(data_depths)
data_depths <- data_depths %>%
  group_by(chr)%>% 
  dplyr::summarise(across(everything(), list(median)))
colnames(data_depths) <- names
chr <- data_depths$chr
data_depths <- as.matrix(data_depths[ , -1])
data_depths_median <- colMedians(data_depths)
data_depths_MAD <- colMads(data_depths)
data_depths <- (data_depths - data_depths_median)/data_depths_MAD
row.names(data_depths) <- chr
data_depths <- data_depths[row.names(data_depths) %in% c("3", "6", "8"), ]

### Create a healthy median
healthy_depths <- data_depths[ , !(colnames(data_depths) %in% data_samples$sWGS)]
data_depths <- data_depths[ , colnames(data_depths) %in% data_samples$sWGS]
healthy_means <- rowMeans2(healthy_depths)
healthy_sd <- rowSds(healthy_depths)

### Convert to Z-scores
data_depths <- (data_depths - healthy_means)/healthy_sd
sample_z <- colSums(abs(data_depths))
data_samples <- merge(data_samples, sample_z, by.x = "sWGS", by.y = "row.names", all = TRUE)

healthy_zscores <- colSums2(abs(healthy_depths))
z_limit <- quantile(healthy_zscores, 0.9)

### Generate Z-score statistics for ichorCNA data (short)
### Remove problem bins and NA bins
data_depths_short <- merge(data_depths_short, healthy_depths_short, by = c("chr", "start", "end"), all = TRUE)
problems <- c(285, 434, 677, 710, 711, 983, 984, 1187, 1263)
data_depths_short <- data_depths_short[!(rownames(data_depths_short) %in% problems), -c(2,3)]
data_depths_short <- data_depths_short[complete.cases(data_depths_short), ]

### Normalize copy number across chromosomes (short)
names <- colnames(data_depths_short)
data_depths_short <- data_depths_short %>%
  group_by(chr)%>% 
  dplyr::summarise(across(everything(), list(median)))
colnames(data_depths_short) <- names
chr <- data_depths_short$chr
data_depths_short <- as.matrix(data_depths_short[ , -1])
data_depths_short_median <- colMedians(data_depths_short)
data_depths_short_MAD <- colMads(data_depths_short)
data_depths_short <- (data_depths_short - data_depths_short_median)/data_depths_short_MAD
row.names(data_depths_short) <- chr
data_depths_short <- data_depths_short[row.names(data_depths_short) %in% c("3", "6", "8"), ]

### Create a healthy median (short)
healthy_depths_short <- data_depths_short[ , !(colnames(data_depths_short) %in% data_samples$sWGS)]
data_depths_short <- data_depths_short[ , colnames(data_depths_short) %in% data_samples$sWGS]
healthy_means <- rowMeans2(healthy_depths_short)
healthy_sd <- rowSds(healthy_depths_short)

### Convert to Z-scores (short)
data_depths_short <- (data_depths_short - healthy_means)/healthy_sd
sample_z_short <- colSums(abs(data_depths_short))
data_samples <- merge(data_samples, sample_z_short, by.x = "sWGS", by.y = "row.names", all = TRUE)

healthy_zscores_short <- colSums2(abs(healthy_depths_short))
z_limit_short <- quantile(healthy_zscores_short, 0.9)

### Rename Z-score columns in metadata
names(data_samples)[names(data_samples) == "y.x"] <- "zscore"
names(data_samples)[names(data_samples) == "y.y"] <- "zscore_short"

### Format annotations and matrices for oncoplot
### Remove normal buffy coats/HBC controls and order samples and oncoplot matrix
data_samples <- data_samples[!(data_samples$Patient == "UMB-HBC" |
                                 data_samples$Timepoint == "Lymphocytes"), ]
data_samples$Patient <- gsub("UMB-0", "", data_samples$Patient)
data_samples <- data_samples[order(factor(data_samples$Relapse, levels = c("Yes", "No", ""),
                                          labels = c("Yes", "No", "Lost to Follow-up")),
                             factor(data_samples$Stage, levels = c("IIB", "IIIA", "IIIB")),
                             factor(data_samples$Treatment, levels = c("Brachytherapy", "Enucleation")),
                             factor(data_samples$Volume),
                             data_samples$Patient,
                             factor(data_samples$Timepoint, levels = c("Tumour", "Baseline", "2 weeks", "3 months", "6 months", "12 months"))), ]
row.names(data_onco) <- data_onco$file_ID
data_onco <- data_onco[data_onco$file_ID %in% data_samples$Targeted.Panel, -c(1:3)]
data_samples <- data_samples[data_samples$Targeted.Panel %in% row.names(data_onco), ]
data_onco <- data_onco[data_samples$Targeted.Panel, c("Chr3", "Chr6p", "Chr6q", "Chr8p", "Chr8q", "GNA11", "GNAQ", "BAP1", "SF3B1", "EIF1AX")]
colnames(data_onco) <- c("Chr 3", "Chr 6p", "Chr 6q", "Chr 8p", "Chr 8q", "GNA11", "GNAQ", "BAP1", "SF3B1", "EIF1AX")

### Set Z-score annotation
data_z <- data_samples[ , c("Targeted.Panel", "zscore", "zscore_short")]
data_z$limit <- z_limit
data_z$limit_short <- z_limit_short
row.names(data_z) <- data_z$Targeted.Panel
data_z <- as.matrix(data_z[ , 2:5])
class(data_z) <- "numeric"
colnames(data_z) <- c("All Fragments", "Short Fragments", "All Limit", "Short Limit")

data_z_short <- data_z[, c("Short Fragments", "Short Limit")]
data_z <- data_z[, c("All Fragments", "All Limit")]

### Set ichorCNA annotation
data_samples$line <- 3
data_samples$ichorCNA <- ifelse(data_samples$TF > data_samples$TF_short, data_samples$TF, data_samples$TF_short)*100
data_ichor <- as.matrix(data_samples[ , c("ichorCNA", "line")])
row.names(data_ichor) <- data_samples$Targeted.Panel

### Set clinical annotations
data_stage <- as.matrix(data_samples$Stage)
row.names(data_stage) <- data_samples$Targeted.Panel

data_treatment <- as.matrix(data_samples$Treatment)
row.names(data_treatment) <- data_samples$Targeted.Panel

data_relapse <- as.matrix(data_samples$Relapse)
row.names(data_relapse) <- data_samples$Targeted.Panel

data_size <- as.matrix(data_samples$Volume)
row.names(data_size) <- data_samples$Targeted.Panel
volume <- expression("Tumour Volume mm"^3)

data_sex <- as.matrix(data_samples$Sex)
row.names(data_sex) <- data_samples$Targeted.Panel

data_samples$Age <- ifelse((data_samples$Age > 0 & data_samples$Age < 50), "0-49",
                        ifelse((data_samples$Age > 49 & data_samples$Age < 60), "50-59",
                               ifelse((data_samples$Age > 59 & data_samples$Age < 70), "60-69", "70-79")))
data_age <- as.matrix(data_samples$Age)
row.names(data_age) <- data_samples$Targeted.Panel
data_age <- factor(data_age, levels = c("0-49", "50-59", "60-69", "70-79"))

data_labels <- as.matrix(data_samples$Timepoint)
row.names(data_labels) <- data_samples$Targeted.Panel

### Calculate alteration percentages for annotation
data_tumours <- data_samples[data_samples$Timepoint == "Tumour", ]
data_pct <- data_onco[row.names(data_onco) %in% data_samples$Targeted.Panel, ]
genes <- colnames(data_pct)
data_pct[data_pct == ""] <- NA 
data_pct[data_pct == "missense"] <- 1
data_pct[data_pct == "frameshift"] <- 1
data_pct[data_pct == "LOH"] <- 1
data_pct[data_pct == "gain"] <- 1
data_pct <- mutate_all(data_pct, function(x) as.numeric(as.character(x)))
data_pct <- trunc((colSums(data_pct, na.rm=TRUE)/11)*100)
data_pct <- paste0(data_pct, '%')
data_pct <- as.matrix(data_pct)
rownames(data_pct) <- genes

### Transpose oncoplot matrix
data_onco <- as.data.frame(t(data_onco))
data_onco <- as.matrix(data_onco)

### Set annotation and oncoplot colours
col <- c(missense = "#33a02c", missense_plasma = "#b2df8a",
         frameshift = "black", frameshift_plasma = "gray30",
         inframe_deletion = "#FF7F00", inframe_deletion_plasma = "#FDBF6F",
         splice_site = "#6A3D9A", splice_site_plasma = "#CAB2D6",
         NS = "black", 
         LOH = "#1f78b4", LOH_plasma = "#a6cee3", 
         gain = "#e31a1c", gain_plasma = "#fb9a99")
col_stage <- c(IIB = "#FFFFB3", IIIA = "#FDB462", IIIB = "#fb9a99")
col_treatment <- c(Brachytherapy = "#FFFFB3", Enucleation = "#FDB462")
col_relapse <- c(Yes = "#fb9a99", No = "#a6cee3", "Lost to Follow-up" = "lightgrey")
col_size <- colorRamp2(c(500, 3000), c("#A1D99B", "#00441B"))
col_age <- c("0-49" = "#A1D99B", "50-59" = "#6BA770", "60-69" = "#00441B", "70-79" = "#00441B")
col_sex <- c(Male = "#a6cee3", Female = "#fb9a99")

### Set variables
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = "grey95", col = NA)),
  missense = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["missense"], col = NA)),
  frameshift = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["frameshift"], col = NA)),
  LOH = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["LOH"], col = NA)),
  gain = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["gain"], col = NA)),
  inframe_deletion = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["inframe_deletion"], col = NA)),
  splice_site = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["splice_site"], col = NA)),
  missense_plasma = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["missense_plasma"], col = NA)),
  frameshift_plasma = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["frameshift_plasma"], col = NA)),
  LOH_plasma = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["LOH_plasma"], col = NA)),
  gain_plasma = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["gain_plasma"], col = NA)),
  inframe_deletion_plasma = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["inframe_deletion_plasma"], col = NA)),
  splice_site_plasma = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["splice_site_plasma"], col = NA)),
  NS = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.1, gp = gpar(fill = col["NS"], col = NA)))

### Compile annotation layers
top_annotation <- HeatmapAnnotation("IchorCNA\nTF (%)" = anno_lines(data_ichor, 
                                                                    add_points = TRUE,
                                                                    pch = c(16, NA),
                                                                    gp = gpar(col = c("black", "red"),
                                                                              lty = c("solid", "dashed")),
                                                                    ylim = c(-1, 10),
                                                                    axis_param = list(gp = gpar(fontsize = 10))),
                                    "Log2(Z-Score)\n(Chr 3, 6, 8)\nAll Fragments" = anno_lines(log2(data_z), 
                                                                                               add_points = TRUE,
                                                                                               pch = c(16, NA),
                                                                                               gp = gpar (col = c("#1f78b4", "#1f78b4"),
                                                                                                          lty = c("solid", "dashed")),
                                                                                               pt_gp = gpar(col = c("#1f78b4", "#1f78b4")),
                                                                                               ylim = c(-1, 6),
                                                                                               axis_param = list(gp = gpar(fontsize = 10))),
                                    "Log2(Z-Score)\n(Chr 3, 6, 8)\nShort Fragments" = anno_lines(log2(data_z_short), 
                                                                                                 add_points = TRUE,
                                                                                                 pch = c(16, NA),
                                                                                                 gp = gpar (col = c("#e31a1c", "#e31a1c"),
                                                                                                            lty = c("solid", "dashed")),
                                                                                                 pt_gp = gpar(col = c("#e31a1c", "#e31a1c")),
                                                                                                 ylim = c(-1, 6),
                                                                                                 axis_param = list(gp = gpar(fontsize = 10))),
                                    "Sex" = data_sex,
                                    "Age" = data_age,
                                    "Tumour Volume" = data_size,
                                    "Treatment" = data_treatment,
                                    "Stage" = data_stage,
                                    "Relapse" = data_relapse,
                                    annotation_name_side = "left",
                                    annotation_name_gp = gpar(fontsize = 13),
                                    height = unit(10, "cm"),
                                    simple_anno_size = unit(0.5, "cm"),
                                    border = TRUE,
                                    annotation_name_rot = 0,
                                    col = list("Age" = col_age, "Sex" = col_sex, "Stage" = col_stage, "Tumour Volume" = col_size, 
                                               "Treatment" = col_treatment, "Relapse" = col_relapse))

label_annotation <- HeatmapAnnotation(Timepoint = anno_text(data_labels,
                                                            gp = gpar(fontsize = 13)))

pct_annotation = rowAnnotation(new_pct = anno_text(data_pct, 
                                                   gp = gpar(fontsize = 13)))

### Set labels and legends
annotation_legend = packLegend(list = list(Legend(title = "CNA Z-Score",
                                type = "points",
                                at = c("All Fragments", "Short Fragments"),
                                labels = c("All Fragments", "Short Fragments"),
                                legend_gp = gpar(col = c("#1f78b4", "#e31a1c")),
                                background = "white"),
                         Legend(title = "Sex", 
                                at = c("Female", "Male"),
                                legend_gp = gpar(fill = c("#fb9a99", "#a6cee3"))),
                         Legend(title = "Age",
                                at = c("0-49", "50-59", "60-69", "70-79"),
                                legend_gp = gpar(fill = c("#A1D99B", "#6BA770", "#00441B", "#00441B"))),
                         Legend(title = expression(bold("Tumour Volume mm"^3)),
                                col_fun = colorRamp2(c(500, 3000), c("#A1D99B", "#00441B"))),
                         Legend(title = "Treatment",
                                at = c("Brachytherapy", "Enucleation"),
                                legend_gp = gpar(fill = c("#FFFFB3", "#FDB462"))),
                         Legend(title = "Tumour Stage",
                                at = c("IIB", "IIIA", "IIIB"),
                                legend_gp = gpar(fill = c("#FFFFB3", "#FDB462", "#fb9a99"))),
                         Legend(title = "Outcome",
                                at = c("Relapsed", "Remission", "Lost to Follow-up"),
                                legend_gp = gpar(fill = c("#fb9a99", "#a6cee3", "lightgrey"))),
                         Legend(title = "Alterations",
                                type = "lines",
                                at = c("Missense Tumour", "Frameshift Tumour", "Inframe Deletion Tumour", "Splice Site Tumour", "LOH Tumour", "Gain Tumour",
                                       "Missense Plasma", "Frameshift Plasma", "LOH Plasma", "Gain Plasma", "Not Sequenced"),
                                legend_gp = gpar(col = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "black"),
                                                 size = 2),
                                background = c("#33a02c", "black", "#FF7F00", "#6A3D9A", "#1f78b4", "#e31a1c", 
                                               "#b2df8a", "gray30", "#a6cee3", "#fb9a99", "grey"))
                         ))

row_labels = row.names(data_onco)

### Set column/row order and splits
sample_order <- data_samples$Targeted.Panel
gene_order <- row.names(data_onco)
column_split <- data_samples$Patient
split_order <- unique(column_split)
column_split <- factor(column_split, levels = split_order)

### Generate oncoprint
pdf(file.path(outdir, "Figure 1.pdf"), width = 12, height = 8)
oncoPrint <- oncoPrint(data_onco,
                       alter_fun = alter_fun, 
                       col = col,
                       row_names_side = "left", pct_side = "none",
                       show_heatmap_legend = FALSE,
                       top_annotation = top_annotation,
                       bottom_annotation = label_annotation,
                       right_annotation = pct_annotation,
                       row_order = gene_order,
                       column_order = sample_order,
                       row_labels = row_labels,
                       row_names_gp = gpar(fontsize = 13),
                       column_names_gp = gpar(fontsize = 13),
                       column_split = column_split,
                       border = TRUE,
                       border_gp = gpar(col = "black"))
draw(oncoPrint, heatmap_legend_list = annotation_legend, show_annotation_legend = FALSE)
dev.off()

### Prepare metrics to export
data_vaf <- file.path(path, "pipeline_analysis", "variants_vaf.txt")
data_vaf <- read.delim(data_vaf)
data_vaf[is.na(data_vaf)] <- 0
data_vaf$vaf <- rowMaxs(as.matrix(data_vaf[, 2:6]))
data_vaf$vaf <- data_vaf$vaf/colMaxs(as.matrix(data_vaf$vaf))
data_vaf <- data_vaf[, colnames(data_vaf) %in% c("sample", "vaf")]
write.table(data_vaf, file.path(outdir, "variant_vafs.txt"), sep = "\t")

data_zscores <- cbind(as.data.frame(data_z[complete.cases(data_z), ]),
                      as.data.frame(data_z_short[complete.cases(data_z_short), ]))
data_zscores$sample <- row.names(data_zscores)
data_zscores$ichor <- rowMaxs(as.matrix(data_zscores[, c(1,3)]))
data_zscores$ichor <- log(data_zscores$ichor)
data_zscores$limit <- ifelse(z_limit > z_limit_short, z_limit, z_limit_short)
data_zscores$limit <- log(data_zscores$limit)

data_zscores$limit <- data_zscores$limit/colMaxs(as.matrix(data_zscores$ichor))
data_zscores$ichor <- abs(data_zscores$ichor)/colMaxs(as.matrix(data_zscores$ichor))
data_zscores <- data_zscores[, c("sample", "ichor", "limit")]
write.table(data_zscores, file.path(outdir, "ichor_score.txt"), sep = "\t", row.names = FALSE)
