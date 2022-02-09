library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project"
healthy_dir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/combined_cohort"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Uveal Melanoma/Manuscript/Figures/Supplemental Figure 5"

data_frequency <- file.path(path, "insert_size", "TGL48_UVM_fragment_swgs_freq.txt")
data_samples <- file.path(path, "TGL48_UVM_sample_list.txt")

healthy_frequency <- file.path(healthy_dir, "insert_size/HBC_fragment_freq.txt")

### Read in data
data_frequency <- read.delim(data_frequency)
data_samples <- read.delim(data_samples)
healthy_frequency <- read.delim(healthy_frequency)

### Format sample sheet
data_samples$Relapse <- factor(data_samples$Relapse, levels = c("No", "Yes", ""),
                               labels = c("No", "Yes", "LtFU"))

### Format dataframes
data_frequency[is.na(data_frequency)] <- 0
data_frequency <- merge(data_frequency, healthy_frequency, by = "length")
row.names(data_frequency) <- data_frequency$length

### Separate uveal and normal samples from frequency distributions
data_uveal <- data_frequency[ , colnames(data_frequency) %in% data_samples$sWGS]
data_uveal$length <- data_frequency$length
data_normal <- data_frequency[ , !(colnames(data_frequency) %in% data_samples$sWGS) ]

### Melt frequency data and append metadata
data_melt <- melt(data_uveal,  id.vars = 'length', variable.name = 'sample')
data_melt$length <- as.numeric(data_melt$length)
data_melt$value <- as.numeric(data_melt$value)
data_melt <- merge(data_melt, data_samples, by.x = c("sample"), by.y = c("sWGS"))

### Make healthy normal frequency distribution
data_normal <- data_normal[, -c(1)]
data_normal <- as.numeric(rowMeans(data_normal))
data_normal <- as.data.frame(data_normal)
data_normal$length <- as.numeric(data_uveal$length)

## Graph Figure 2B - Fragment Frequencies
Fig <- ggplot(data_melt, aes(length, value)) +
  geom_line(aes(group = sample), alpha = 0.3) +
  geom_line(data = data_normal, aes(length, data_normal, color = "#fb9a99"), size = 0.75, alpha = 0.75) +
  facet_wrap(.~Relapse) +
  xlab("") + 
  ylab("% Frequency") +
  ggtitle("Relapse Status") + 
  labs(color = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 15), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 15),
        legend.position = "none",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 1)) + 
  scale_y_continuous(limits=c(0, 3.5), expand = c(0,0)) + 
  scale_x_continuous(limits=c(0,320, expand = c(0,0))) +
  geom_vline(xintercept=167, linetype="dashed", color = "black", size=0.75)
Fig

subFig <- ggplot(data_melt, aes(length, value)) +
  geom_line(aes(group = sample), alpha = 0.3) +
  geom_line(data = data_normal, aes(length, data_normal, color = "#fb9a99"), size = 0.75, alpha = 0.75) +
  facet_wrap(.~Relapse) +
  xlab("Fragment Size (bp)") + 
  ylab("% Frequency") +
  ggtitle("") + 
  labs(color = "") +
  scale_color_discrete(labels = c("Healthy")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 15),
        legend.position = "none",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 1)) + 
  scale_y_continuous(limits=c(0, 0.5), expand = c(0,0)) + 
  scale_x_continuous(limits=c(40,150, expand = c(0,0))) +
  geom_vline(xintercept=167, linetype="dashed", color = "black", size=0.75)
subFig

gA <- ggplotGrob(Fig)
gB <- ggplotGrob(subFig)

pdf(file.path(outdir, "Supplemental Figure 5.pdf"), width = 8, height = 6)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))
dev.off()
