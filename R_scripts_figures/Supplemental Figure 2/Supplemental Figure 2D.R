library(tidyverse)
library(reshape2)
library(matrixStats)
library(GenomicRanges)
library(grid)

### Set variables
data_ratio <- "TGL48_UVM_fragment_ratio_5Mb.txt"
samples <- "TGL48_UVM_sample_list.txt"
outdir <- ""

### Read in data
data_ratio <- read.delim(data_ratio, check.names = FALSE)
samples <- read.delim(samples, check.names = FALSE)

### Subset plasma samples from uveal data
samples <- samples[!(samples$Timepoint %in% c("Lymphocytes", "Tumour", "Healthy")), ]


### Seperate out Chr15
data_ratio <- data_ratio[data_ratio$seqnames == "chr15", ]
data_ratio <- data_ratio[, samples$sWGS]

### Melt data
data_melt <- melt(data_ratio)
data_melt <- merge(data_melt, samples, by.x = "variable", by.y = "sWGS", all = TRUE)

### Format melted data
data_melt$Patient <- sub("UMB-0", "", data_melt$Patient)
order <- unique(data_melt$Patient)
data_melt$Patient <- factor(data_melt$Patient, levels = order)

### Plot data
ratio_plot <- ggplot(data_melt) + 
  geom_boxplot(aes(Patient, value, fill = Patient), alpha = 0.5) +
  ggtitle("Chr15") +
  xlab("Patient") +
  ylab("Fragment Ratio") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 10),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 10)) + 
  scale_fill_manual(values = c(rep("grey", 9), "#fa4e0c", "grey")) +
  scale_y_continuous(limits = c(-0.04, 0.04))
ratio_plot
ggsave(file.path(outdir, paste0("chr15_ratios.pdf")), ratio_plot, device = "pdf", width = 3.5, height = 3, units = "in")
