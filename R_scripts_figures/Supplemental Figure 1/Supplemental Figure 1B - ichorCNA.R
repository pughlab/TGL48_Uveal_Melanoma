library(ggplot2)
library(reshape2)
library(plyr)
library(ggforce)

### Set variables
path <- ""
outdir <- ""

data_samples <- file.path(path, "TGL48_UVM_sample_list.txt")

### Import data
data_samples <- read.delim(data_samples)

### Keep plasma samples only and format dataframe
data_samples <- data_samples[!(data_samples$Timepoint %in% c("Tumour", "Lymphocytes", "Healthy")), ]
data_samples$Timepoint <- as.character(data_samples$Timepoint)
data_samples$Timepoint <- factor(data_samples$Timepoint, 
                                 levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))
data_samples$Patient <- gsub("UMB-0", "", data_samples$Patient)
order <- unique(data_samples$Patient)
data_samples$Patient <- factor(data_samples$Patient, levels = order)

### Remove short TFs with coverage > 0.1x and melt data
data_samples$TF_short <- ifelse(data_samples$coverage_short >= 0.1, data_samples$TF_short, NA)
data_samples <- data_samples[, c("Patient", "Timepoint", "TF", "TF_short")]
data_samples <- melt(data_samples)
data_samples$variable <- factor(data_samples$variable, 
                               levels = c("TF", "TF_short"), 
                               labels = c("all", "90-150bp"))
data_samples$detection <- ifelse(data_samples$value >= 0.03, "yes", "no")

### Plot Supplemental Figure 1B
tumour_fraction <- ggplot(data_samples, aes(x = Timepoint, 
                                            y = value,
                                            group = Patient,
                                            fill = detection)) +
  geom_point(aes(), alpha = 0.8, size = 2, shape = 21) +
  xlab("Timepoint") + 
  ylab("ichorCNA Tumour Fraction") +
  scale_fill_manual(values = c("grey", "#CC2D35")) +
  facet_grid(variable ~ Patient, 
             scales = "free", 
             space = "free") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        strip.background = element_rect(color="white", fill="white"),
        strip.text.x = element_text(size = 12, color = "black", face = "bold", angle = 0),
        strip.text.y = element_text(size = 12, color = "black", face = "bold", angle = 270)) + 
  scale_y_continuous(limits=c(-0.01, 0.1), expand = c(0,0)) +
  geom_hline(data = data_samples, aes(yintercept = 0.03), linetype = "dashed", color = "grey")
tumour_fraction

ggsave(file.path(outdir, "Supplemental Figure 3 - ichorCNA.pdf"), tumour_fraction, device = "pdf", width = 8, height = 5, units = "in")
