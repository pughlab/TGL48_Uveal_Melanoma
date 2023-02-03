library(tidyverse)
library(dplyr)
library(data.table)
library(reshape2)
library(lemon)
library(matrixStats)

### Set paths
path <- ""
outdir <- ""
PRC1 <- "H2AK119Ub_targets.txt"
healthy <- "HBC_fragment_ratios_100kb.txt"

frags <- file.path(path, "/fragmentomics/TGL48_UVM_fragment_ratio_100kb.txt")
samples <- file.path(path, "TGL48_UVM_sample_list.txt")

### Import data
data_frag <- read.delim(frags)
data_prc1 <- read.delim(PRC1)
samples <- read.delim(samples)
HBC <- read.delim(healthy, check.names = FALSE)

### Format sample list
samples <- samples[!(samples$Timepoint %in% c("Lymphocytes", "Tumour", "Healthy")), ]
samples$Patient <- gsub("UMB-0", "", samples$Patient)
samples <- samples[samples$Timepoint %in% c("Baseline", "2 weeks", "3 months", "6 months", "12 months"), ]
samples$Relapse <- factor(samples$Relapse, levels = c("No", "Yes", ""),
                          labels = c("No", "Yes", "Lost to Follow-up"))
samples$timepoint <- as.character(samples$Timepoint)
samples$timepoint <- factor(samples$Timepoint,
                            levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))
samples <- samples[order(factor(samples$Relapse, levels = c("Yes", "No", "Lost to Follow-up")),
                         factor(samples$Stage, levels = c("IIB", "IIIA", "IIIB")),
                         factor(samples$Treatment, levels = c("Brachytherapy", "Enucleation")),
                         samples$Patient,
                         factor(samples$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))), ]

### Normalize and center 100kb bin ratios (Uveal cohort)
chr <- data_frag[ , c(1:3)]
data_frag <- as.matrix(data_frag[ , -c(1:3)])
data_frag <- ((data_frag - colMeans2(data_frag, na.rm = TRUE))/colSds(data_frag, na.rm = TRUE))*0.01
data_frag <- cbind(chr, data_frag)

### Normalize and center 100kb bin ratios (HBC cohort)
HBC <- as.matrix(HBC[ , -c(1:3)])
HBC <- ((HBC - colMeans2(HBC, na.rm = TRUE))/colSds(HBC, na.rm = TRUE))*0.01
HBC <- cbind(chr, HBC)

### Find overlapping bins
colnames(data_prc1) <- c("ensemble", "gene", "seqnames", "start", "end")
prc1_targets <- as.data.table(data_prc1[ , c(3,4,5)])
frag_targets <- as.data.table(chr)

setDT(prc1_targets)
setkey(frag_targets)
overlap <- foverlaps(prc1_targets, frag_targets, type="within", nomatch=0L)

### Merge overlap with ratios
overlap <- overlap[ , c(1:3)]
overlap <- overlap[!duplicated(overlap),]
overlap$PRC1 <- "yes"
data_frag <- merge(data_frag, overlap, by = c("seqnames", "start", "end"), all = TRUE)
data_frag$PRC1[is.na(data_frag$PRC1)] <- "no"

HBC <- merge(HBC, overlap, by = c("seqnames", "start", "end"), all = TRUE)
HBC$PRC1[is.na(HBC$PRC1)] <- "no"

### Format ratio dataframe
data_frag <- data_frag[data_frag$PRC1 == "yes", ]
data_frag <- data_frag[ , samples$sWGS]
data_frag <- data_frag[complete.cases(data_frag), ]

### Calculate Healthy median and apphend to data
HBC <- as.matrix(HBC[HBC$PRC1 == "yes", !(colnames(HBC) %in% c("seqnames", "start", "end", "PRC1"))])
HBC <- HBC[row.names(HBC) %in% row.names(data_frag), ]
HBC_median <- rowMedians(HBC, na.rm = TRUE)

data_frag$HBC <- HBC_median
data_frag <- as.matrix(data_frag)
HBC_sample <- c(replicate(5, "HBC"), NA, NA, replicate(5, "HBC"))
samples <- rbind(samples, HBC_sample)

### Calculate Mean and SD of PRC1 target loci
mean <- colMeans2(data_frag)
sd <- colSds(data_frag)
summary <- cbind(mean, sd)
row.names(summary) <- colnames(data_frag)
summary <- summary[!(rownames(summary) == "HBC"), ]
summary <- merge(summary, samples, by.x = "row.names", by.y = "sWGS")
summary$Relapse <- factor(summary$Relapse, levels = c("Yes", "No", "Lost to Follow-up"))
summary$Timepoint <- factor(summary$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))

### Melt data
data_melt <- melt(data_frag)
data_melt <- merge(data_melt, samples, by.x = "Var2", by.y = "sWGS")
order <- unique(samples$Patient)
data_melt$Patient <- factor(data_melt$Patient, levels = order)
data_melt$Relapse[is.na(data_melt$Relapse)] <- "HBC"
data_melt$Relapse <- factor(data_melt$Relapse, levels = c("HBC", "Yes", "No", "Lost to Follow-up"))
data_melt$Timepoint <- factor(data_melt$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))

### Plot Summary
plot <- ggplot(summary) + 
  geom_point(aes(mean, sd, color = Relapse)) +
  labs(color="Relapse") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 10),
        legend.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) + 
  scale_color_manual(values = c("#fa4e0c", "#3f9e80", "#516eb0", "black"))
plot
ggsave(file.path(outdir, paste0("Supplemental Figure 8 - Means.pdf")), plot, device = "pdf", width = 4.5, height = 3.5, units = "in")

### Plot Ratios
ratio_plot <- ggplot(data_melt) + 
  geom_boxplot(aes(Timepoint, value, fill = Relapse), alpha = 0.5, outlier.shape = NA) +
  facet_grid(.~Patient, scales = "free", space = "free") +
  labs(color="Relapse") +
  xlab("Timepoint") +
  ylab("Fragment Ratio") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
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
        legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 10)) + 
  scale_fill_manual(values = c("grey", "#fa4e0c", "#3f9e80", "#516eb0")) +
  scale_y_continuous(limits = c(0, 0.075))
ratio_plot
ggsave(file.path(outdir, paste0("Supplemental Figure 9 - Ratios.pdf")), ratio_plot, device = "pdf", width = 8, height = 3.5, units = "in")

### Compare PRC1 ratios to genome wide ratios
ratios <- read.delim("/Volumes/GoogleDrive/My Drive/Post-Doc/Uveal Melanoma/Manuscript/Figures/Figure 3/ratios.txt")
data_scores <- read.delim("/Volumes/GoogleDrive/My Drive/Post-Doc/Uveal Melanoma/Manuscript/Figures/Figure 3/PRC1_scores.txt")
ratios <- merge(ratios, data_scores, by = "sample")
ratios <- merge(ratios, samples, by.x = "sample", by.y = "sWGS")
ratios$Relapse <- factor(ratios$Relapse, levels = c("No", "Lost to Follow-up", "Yes"))
ratios <- ratios[order(factor(ratios$Relapse)), ]
order <- unique(ratios$Patient)
ratios$Patient <- factor(ratios$Patient, levels = order)

scores_plot <- ggplot(ratios) + 
  geom_point(aes(ratio, PRC1, color = Relapse)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Genome Wide Score") +
  ylab("PRC1 Target Loci Score") +
  labs(color="Relapse") +
  facet_wrap(.~Patient, nrow = 2) +
  facet_rep_wrap(.~Patient, nrow = 2) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
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
        legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 10)) + 
  scale_color_manual(values = c("#3f9e80", "#516eb0", "#fa4e0c"))
scores_plot
ggsave(file.path(outdir, paste0("Supplemental Figure 10 - Scores.pdf")), scores_plot, device = "pdf", width = 8, height = 4, units = "in")
