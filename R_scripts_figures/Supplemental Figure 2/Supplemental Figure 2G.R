library(tidyverse)
library(dplyr)
library(data.table)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)
library(cowplot)
library(vroom)

### Set paths
path <- ""
outdir <- ""
PRC1 <- "H2AK119Ub_targets.txt"
healthy <- "HBC_fragment_ratios_100kb.txt"
healthy_samples <- "HBC_sample_list.txt"

frags <- file.path(path, "/fragmentomics/TGL48_UVM_fragment_ratio_100kb.txt")
samples <- file.path(path, "TGL48_UVM_sample_list.txt")

### Import data
data_frag <- vroom(frags)
data_prc1 <- read.delim(PRC1)
samples <- read.delim(samples)
HBC <- vroom(healthy)
HBC_samples <- read.delim(healthy_samples)

### Format HBCs
#HBC_samples <- HBC_samples[!(HBC_samples$cohort %in% c("Phallen")), ]
#HBC <- HBC[, colnames(HBC) %in% c("seqnames", "start", "end", HBC_samples$sample)]

### Format sample list
samples <- samples[!(samples$Timepoint %in% c("Lymphocytes", "Tumour", "Healthy")), ]
samples$Patient <- gsub("UMB-0", "", samples$Patient)
samples <- samples[samples$Timepoint %in% c("Baseline", "2 weeks", "3 months", "6 months", "12 months"), ]
samples$Relapse <- factor(samples$Relapse, levels = c("No", "Yes", ""),
                          labels = c("No", "Yes", "LtFU"))
samples$timepoint <- as.character(samples$Timepoint)
samples$timepoint <- factor(samples$Timepoint,
                            levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))
samples <- samples[order(factor(samples$Patient, levels = c("09", "03", "05", "04", "06", "01", "02", "11", "07", "08", "10")),
                         factor(samples$Relapse, levels = c("Yes", "No", "LOF")),
                         factor(samples$Stage, levels = c("IIB", "IIIA", "IIIB")),
                         factor(samples$Treatment, levels = c("Brachytherapy", "Enucleation")),
                         samples$Patient,
                         factor(samples$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))), ]

### Normalize and center 100kb bin ratios (Uveal cohort)
chr <- data_frag[ , c(1:3)]
data_frag <- as.matrix(data_frag[ , -c(1:3)])
data_frag[data_frag > 5] <- 5
data_frag[data_frag < -5] <- -5
data_frag <- scale(data_frag, center = TRUE)
data_frag <- cbind(chr, data_frag)

### Normalize and center 100kb bin ratios (HBC cohort)
HBC <- as.matrix(HBC[ , -c(1:3)])
HBC <- scale(HBC, center = TRUE)
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

### Format ratio dataframe (select non PRC1 target bins)
data_frag <- data_frag[data_frag$PRC1 == "no", ]
data_frag <- data_frag[ , samples$sWGS]
data_frag <- data_frag[complete.cases(data_frag), ]

### Make sure sample order is same as sample sheet
data_frag <- data_frag[ , colnames(data_frag) %in% samples$sWGS]
data_frag <- data_frag[ , samples$sWGS]

### Find healthy median of PRC1 Target Loci
HBC <- as.matrix(HBC[HBC$PRC1 == "no", !(colnames(HBC) %in% c("seqnames", "start", "end", "PRC1"))])
HBC <- HBC[row.names(HBC) %in% row.names(data_frag), ]
HBC_median <- rowMedians(HBC, na.rm = TRUE)
HBC_sds <- rowSds(HBC, na.rm = TRUE)

### Select random bins (same # as PRC1 bins)
data_frag <- data_frag[sample(nrow(data_frag), nrow(overlap)), ]
HBC <- HBC[row.names(data_frag), ]

### Calculate HBC Z-scores
mean <- mean(HBC_median)
sd <- sd(HBC_median)
HBC_means <- colMeans2(as.matrix(HBC), na.rm = TRUE)
HBC_sds <- colSds(as.matrix(HBC), na.rm = TRUE)
HBC_scores <- (HBC_means - mean)/sqrt((sd^2) + (HBC_sds)^2)
HBC_scores_mean <- mean(HBC_scores)
limit <- quantile(HBC_scores, 0.9)

### Calculate sample Z-scores
data_means <- colMeans2(as.matrix(data_frag), na.rm = TRUE)
data_sds <- colSds(as.matrix(data_frag), na.rm = TRUE)
data_scores <- (data_means - mean)/sqrt((sd^2) + (data_sds)^2)

### Set Z-score annotation
data_scores <- as.data.frame(data_scores)
data_scores$limit <- limit
data_scores <- as.matrix(data_scores)
row.names(data_scores) <- colnames(data_frag)

### Set clinical information
data_scores <- merge(data_scores, samples, by.x = "row.names", by.y = "sWGS")
data_scores$Patient <- factor(data_scores$Patient, levels = c("09", "03", "05", "04", "06", "01", "02", "11", "07", "08", "10"))
data_scores$Timepoint <- factor(data_scores$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))

### Plot
plot <- ggplot(data_scores, aes(Timepoint, data_scores, fill = Relapse)) +
  geom_point(pch = 21) + 
  geom_line(aes(group = Patient)) +
  geom_hline(data = data_scores, aes(yintercept = limit), linetype = "dashed", color = "red") +
  facet_wrap(.~Patient, nrow = 1) +
  scale_fill_manual(values = c( "#3f9e80", "#fa4e0c", "#516eb0")) +
  ggtitle("") +
  xlab("Timepoint") + 
  ylab("Z-score (random bins)") +
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
  scale_y_continuous(expand = c(0,0), limits = c(-1, 2))
plot

ggsave(file.path(outdir, "other_ratio_scores.pdf"), plot, device = "pdf", width = 9, height = 3, units = "in")
