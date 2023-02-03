library(ComplexHeatmap)
library(dplyr)

setwd("")

data_onco <- read.table("Supplemental Figure 1A.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## Set ichorCNA
data_onco$ichorCNA <- data_onco$ichorCNA*100
data_ichor <- data_onco$ichorCNA

## Set alteration percentages
data_pct <- data_onco[, !(colnames(data_onco) %in% c("sample_ID", "detection","ichorCNA"))]
data_pct <- data_pct[c("Chr3", "Chr6p", "Chr6q", "Chr8p", "Chr8q")]
genes <- colnames(data_pct)
data_pct[data_pct == ""] <- NA 
data_pct[data_pct == "LOH"] <- 1
data_pct[data_pct == "gain"] <- 1
data_pct <- mutate_all(data_pct, function(x) as.numeric(as.character(x)))
data_pct <- trunc((colSums(data_pct, na.rm=TRUE)/11)*100)
data_pct <- paste0(data_pct, '%')
data_pct <- as.matrix(data_pct)
rownames(data_pct) <- genes

## Set labels
data_onco$detection <- factor(data_onco$detection, 
                            levels = c("IMPACT", "ichorCNA"), 
                            labels = c("MLPA", "sWGS"))
data_labels <- as.matrix(data_onco$detection)
row.names(data_labels) <- data_onco$sample_ID

## Set samples
data_samples <- data_onco[, c("patient_ID", "sample_ID")]

## Format heatmap
row.names(data_onco) <- data_onco$sample_ID
data_onco <- as.data.frame(t(data_onco))
data_onco <- data_onco[c("Chr3", "Chr6p", "Chr6q", "Chr8p", "Chr8q"), ]
data_onco <- as.matrix(data_onco)

## Set colours
col <- c(LOH = "#1f78b4", LOH_ichor = "#a6cee3", 
         gain = "#e31a1c", gain_ichor = "#fb9a99")

## Set variables
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = "lightgrey", col = NA)),
  LOH = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["LOH"], col = NA)),
  gain = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["gain"], col = NA)),
  LOH_ichor = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["LOH_plasma"], col = NA)),
  gain_ichor = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["gain_plasma"], col = NA)))

## Set additional annotations
top_annotation <- HeatmapAnnotation("sWGS Tumour Fraction (%)" = anno_points(data_ichor, size = unit(3, "mm"), ylim = c(0, 15)),
                                    annotation_name_side = "left",
                                    height = unit(3, "cm"),
                                    annotation_name_rot = 0)

label_annotation <- HeatmapAnnotation(Timepoint = anno_text(data_labels))

freq_annotation <- rowAnnotation(IMPACT_barplot = anno_oncoprint_barplot(c("LOH", "LOH_ichor", "gain", "gain_ichor"),
                                                                         border = TRUE, height = unit(4, "cm"), 
                                                                         axis_param = list(side = "bottom", labels_rot = 0, direction = "reverse")))
                                 #ichor_barplot = anno_oncoprint_barplot(c("LOH_ichor", "gain_ichor"),
                                  #                                      border = TRUE, height = unit(4, "cm"), 
                                   #                                     axis_param = list(side = "bottom", labels_rot = 0, direction = "reverse")))

pct_annotation = rowAnnotation(new_pct = anno_text(data_pct))

## Set labels
heatmap_legend_param = list(title = "Alternations", 
                            at = c("LOH", "gain", "LOH_ichor", "gain_ichor"), 
                            labels = c("LOH MLPA", "Gain MLPA", "LOH sWGS", "Gain sWGS"))

row_labels = structure(c("Chr 3", "Chr 6p", "Chr 6q", "Chr 8p", "Chr 8q"), 
                       names = c("Chr3", "Chr6p", "Chr6q", "Chr8p", "Chr8q"))

## Set order
sample_order <- data_samples$sample_ID
gene_order <- c("Chr3", "Chr6p", "Chr6q", "Chr8p", "Chr8q")
column_split <- data_samples$patient_ID

## Generate oncoprint
oncoPrint <- oncoPrint(data_onco,
                       alter_fun = alter_fun, 
                       col = col,
                       row_names_side = "left", pct_side = "none",
                       top_annotation = top_annotation,
                       bottom_annotation = label_annotation,
                       left_annotation = freq_annotation,
                       right_annotation = pct_annotation,
                       heatmap_legend_param = heatmap_legend_param,
                       row_order = gene_order,
                       column_order = sample_order,
                       row_labels = row_labels,
                       column_split = column_split,
                       border = TRUE,
                       border_gp = gpar(col = "black")
                       )

draw(oncoPrint, heatmap_legend_side = "right")

