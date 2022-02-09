library(tidyverse)
library(reshape2)
library(matrixStats)
library(GenomicRanges)
library(grid)

### Set variables
data_ratio <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project/fragmentomics/TGL48_UVM_fragment_ratio_5Mb.txt"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project/TGL48_UVM_sample_list.txt"
healthy_ratio <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/combined_cohort/fragmentomics/HBC_fragment_ratios.txt"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Uveal Melanoma/Manuscript/Figures/Supplemental Figure 6"

### Read in data
data_ratio <- read.delim(data_ratio, check.names = FALSE)
healthy_ratio <- read.delim(healthy_ratio, check.names = FALSE)
samples <- read.delim(samples, check.names = FALSE)

### Subset plasma samples from uveal data
samples <- samples[!(samples$Timepoint %in% c("Lymphocytes", "Tumour", "Healthy")), ]

### Count number of healthy controls
data_n <- ncol(healthy_ratio) - 4

### Generate healthy median
chromosomes <- healthy_ratio[ , c(1:4)]
healthy <- as.matrix(healthy_ratio[ , -c(1:4)])
healthy_median <- rowMedians(healthy)
healthy_sd <- rowSds(healthy)
healthy_median <- cbind(chromosomes, healthy_median)
healthy_median$lower <- healthy_median$healthy_median - healthy_sd
healthy_median$upper <- healthy_median$healthy_median + healthy_sd
healthy_median$bin <- c(1:nrow(healthy_median))

### Order ratio data and merge together
healthy_ratio$bin <- c(1:nrow(healthy_ratio))
data_ratio$bin <- c(1:nrow(data_ratio))

data_merge <- merge(healthy_ratio, data_ratio, by = c("seqnames", "arm", "start", "end", "bin"))

### Melt data
data_melt <- melt(data_merge, id.vars = c("seqnames", "arm", "start", "end", "bin"))
data_melt <- merge(data_melt, samples, by.x = "variable", by.y = "sWGS", all = TRUE)

### Format melted data
data_melt$Relapse <- ifelse(is.na(data_melt$Relapse), "Healthy", data_melt$Relapse)
data_melt$Relapse <- factor(data_melt$Relapse, levels = c("Healthy", "No", "Yes", ""),
                            labels = c("Healthy", "Remission", "Relapse", "Lost to Follow-up"))
data_melt <- data_melt[order(data_melt$variable,
                             data_melt$bin), ]

# Plot profiles
mytheme <- theme_classic(base_size=12) + theme(
  axis.text.x = element_blank(),
  axis.ticks.x=element_blank(),
  strip.text.x = element_text(size=11),
  strip.text.y = element_text(size=12),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  axis.text.y = element_text(size=15),
  plot.title = element_text(size=15),
  legend.position = "none",
  legend.title = element_text(size=10),
  legend.text = element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_rect(fill="white", color="white"),
  panel.spacing.x=unit(0.1, "lines"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
data_melt$arm <- factor(data_ratio$arm, levels=armlevels)
healthy_median$arm <- factor(healthy_median$arm, levels=armlevels)

arm <- data_ratio %>% group_by(arm) %>%
  dplyr::summarize(n=n()) %>%
  mutate(arm = as.character(arm))
small.arms <- setNames(c("", "10q", "", "12q", "", "16q",
                         "", "17q", "", "18q",
                         "", "", "", "",
                         "", ""),
                       c("10p", "10q", "12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

# Generate Fragmentation and Coverage plots
g1 <- ggplot() + 
  geom_line(data = data_melt, aes(x = bin, y = value, group = variable, color = "red"), size = 0.5, alpha = 0.5) + 
  geom_line(data = healthy_median, aes(x = bin, y = healthy_median), size = 0.75, alpha = 0.5, color = "black") + 
  ggtitle(paste0("")) +
  labs(x="Chromosome", y="Fragmentation profile\n", color="") + 
  facet_grid(Relapse~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels)) + 
  coord_cartesian(xlim = NULL, ylim=c(-.08,.08), expand = TRUE) + 
  mytheme
g1

g2 <- ggplot() + 
  geom_line(data = healthy_median, aes(x = bin, y = healthy_median), size = 0.75, alpha = 0.5, color = "black") + 
  geom_ribbon(data = healthy_median, aes(x = bin, ymin = lower, ymax = upper), fill = "red", alpha = 0.5) + 
  ggtitle(paste0("Healthy Median +/- SD, n= ", data_n)) +
  labs(x="Chromosome", y="Fragmentation profile\n", color="") + 
  facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels)) + 
  coord_cartesian(xlim = NULL, ylim=c(-.08,.08), expand = TRUE) + 
  mytheme
g2

Figure <- ggarrange(g1, g2, ncol = 1, nrow = 2, heights = c(4, 1.5))

ggsave(file.path(outdir, "Supplemental Figure 6.pdf"), Figure, width = 15, height = 9, units = "in")

