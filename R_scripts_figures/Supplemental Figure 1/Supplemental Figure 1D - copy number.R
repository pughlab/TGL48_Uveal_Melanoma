library(ComplexHeatmap)
library(dplyr)
library(reshape2)

path <- ""
outdir <- ""

data_onco <- read.delim(file.path(path, "pipeline_analysis", "tumour_CNV.txt"))

### Perform McNemar's Test
data_test <- data_onco[!(data_onco$Chr3 == "NS"), ]
data_test <- data_test[data_test$Patient %in% data_test$Patient[duplicated(data_test$Patient)],]

data_ichor <- data_test[data_test$detection == "ichorCNA", ]
data_ichor <- reshape2::melt(data_ichor, id = c("Patient", "detection"))
data_impact <- data_test[data_test$detection == "IMPACT", ]
data_impact <- reshape2::melt(data_impact, id = c("Patient", "detection"))

data_test <- merge(data_ichor, data_impact, by = c("Patient", "variable"))
data_test[data_test == ""] <- "FALSE"
data_test[data_test == "LOH" | data_test == "gain"] <- "TRUE"
data_test <- data_test[, colnames(data_test) %in% c("value.x", "value.y")]
colnames(data_test) <- c("ichorCNA", "IMPACT")
test_table <- table(data_test)

mcnemar.test(test_table, y = NULL, correct = TRUE)

## Set labels
data_onco$Patient <- gsub("UMB-0", "", data_onco$Patient)
data_onco$detection <- factor(data_onco$detection, 
                            levels = c("IMPACT", "ichorCNA"))
data_labels <- as.matrix(data_onco$detection)
row.names(data_labels) <- data_onco$sample_ID
data_patient <- data_onco$Patient

## Format heatmap
row.names(data_onco) <- data_onco$sample_ID
data_onco <- as.data.frame(t(data_onco))
data_onco <- data_onco[c("Chr3", "Chr6p", "Chr6q", "Chr8p", "Chr8q"), ]
data_onco <- as.matrix(data_onco)

## Set colours
col <- c(LOH = "#1f78b4", gain = "#e31a1c", NS = "black")
col_method <- c(IMPACT = "grey", ichorCNA = "#B2DF8A")

## Set variables
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = "lightgrey", col = NA)),
  LOH = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["LOH"], col = NA)),
  gain = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["gain"], col = NA)),
  NS = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.1, gp = gpar(fill = col["NS"], col = NA)))

## Set additional annotations
label_annotation <- HeatmapAnnotation(Method = data_labels,
                                      col = list(Method = col_method))

## Set labels
heatmap_legend_param = list(title = "Alternations", 
                            at = c("LOH", "gain", "NS"), 
                            labels = c("LOH", "Gain", "Not Sequenced"))

row_labels = structure(c("Chr 3", "Chr 6p", "Chr 6q", "Chr 8p", "Chr 8q"), 
                       names = c("Chr3", "Chr6p", "Chr6q", "Chr8p", "Chr8q"))

## Set order
sample_order <- colnames(data_onco)
gene_order <- c("Chr3", "Chr6p", "Chr6q", "Chr8p", "Chr8q")
column_split <- data_patient

## Generate oncoprint
pdf(file.path(outdir, "Supplemental Figure 2.pdf"), width = 6, height = 2)
oncoPrint <- oncoPrint(data_onco,
                       alter_fun = alter_fun, 
                       col = col,
                       row_names_side = "left", pct_side = "none",
                       top_annotation = label_annotation,
                       bottom_annotation = NULL,
                       left_annotation = NULL,
                       right_annotation = NULL,
                       heatmap_legend_param = heatmap_legend_param,
                       row_order = gene_order,
                       column_order = sample_order,
                       row_labels = row_labels,
                       column_split = column_split,
                       border = TRUE,
                       border_gp = gpar(col = "black")
                       )

draw(oncoPrint, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()
