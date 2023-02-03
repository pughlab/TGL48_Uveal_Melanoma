library(tidyverse)
library(dplyr)
library(reshape2)
library(ggformula)
library(signal)
library(matrixStats)
library(lemon)

### Set paths
path <- ""
outdir <- ""
healthy <- ""

liver <- file.path(path, "griffin/uveal/TGL48_UVM_griffin_corrected_Liver_dnase.txt")
samples <- file.path(path, "sample_list.txt")

healthy_liver <- file.path(healthy, "HBC_griffin_corrected_Liver_dnase.txt")

### Import data
liver <- read.delim (liver)
samples <- read.delim(samples)
healthy_liver <- read.delim(healthy_liver, check.names = FALSE)

### Format data
distance <- liver$distance
liver <- liver[ , -1]
healthy_liver <- healthy_liver[ , -1]

### Scale data
liver <- as.data.frame(scale(liver)*0.1 + 1)
healthy_liver <- as.data.frame(scale(healthy_liver)*0.1 + 1)

### Subset samples to compare (09 vs 11)
samples <- samples[samples$Patient %in% c("UMB-009", "UMB-011") &
                     !(samples$Timepoint %in% c("Lymphocytes", "Tumour")), ]
liver <- liver[, samples$sWGS]
liver$distance <- distance

### Make healthy median and sd
liver_median <- rowMedians(as.matrix(healthy_liver))
liver_sd <- rowSds(as.matrix(healthy_liver))
liver_median <- as.data.frame(cbind(liver_median, liver_sd))
colnames(liver_median) <- c("median", "sd")
liver_median$distance <- distance
liver_median$variable <- "Healthy"
liver_median$Patient <- "Healthy"
liver_median$type <- "Liver"
liver_median$Timepoint <- "Baseline"

median <- liver_median

### Melt data
liver <- reshape2::melt(liver, id.vars = c("distance"))
liver$type <- "Liver"

### Format melted data
data_melt <- liver
data_melt <- merge(data_melt, samples, by.x = c("variable"), by.y = c("sWGS"))
data_melt$Timepoint <- factor(data_melt$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"),
                              labels = c("Baseline", "2 Weeks", "3 Months", "6 Months", "12 Months"))
data_melt$Patient <- factor(data_melt$Patient, levels = c("UMB-009", "UMB-011"),
                            labels = c("09 - Relapsed", "11 - Remission"))
data_melt <- data_melt[order(data_melt$Patient,
                             data_melt$type,
                             data_melt$Timepoint,
                             data_melt$distance), ]

### Make reference line
data_ref <- data_melt[data_melt$Timepoint == "Baseline", c(2, 4, 5)]
data_ref$value <- 0.8
healthy_ref <- data_ref
healthy_ref$Patient <- "Healthy"
data_ref <- rbind(data_ref, healthy_ref)

### Graph Figure - cohort wide
griffin <- ggplot(data_melt, aes(distance, value)) +
  geom_line(aes(group = variable, color = Timepoint), alpha = 0.75, size = 0.55) +
  geom_line(data = median, aes(distance, median, group = variable, color = Timepoint), alpha = 0.75, size = 0.55) +
  geom_ribbon(data = median, aes(distance, y = median, ymin = (median - sd), ymax = (median + sd), group = variable), alpha = 0.5) + 
  geom_hline(data = data_ref, aes(yintercept = value), linetype = "dashed", alpha = 0.75) +
  scale_color_manual(values = c("Baseline" = "black", "2 Weeks" = "#EF3B2C", "3 Months" = "#CB181D", "6 Months" = "#A50F15", "12 Months" = "#67000D")) + 
  facet_grid(Patient~., scales = "free") +
  facet_rep_grid(Patient~.) + 
  labs(color = "") +
  xlab("Distance") + 
  ylab("Coverage") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.position = "bottom",
        legend.spacing.y = unit(0, "mm"),
        legend.key.size = unit(3, "mm"),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(face = "bold", size = 10),
        panel.spacing = unit(-3, "mm"),
        panel.spacing.y = unit(-5, "mm")) +
  guides(color=guide_legend(nrow = 3,byrow = TRUE)) + 
  scale_x_continuous(limits=c(-500, 500), expand = c(0,0)) +
  scale_y_continuous(limits=c(0.7, 1.2), expand = c(0,0))
griffin
ggsave(file.path(outdir, paste0("Figure 4 - Griffin curves.pdf")), griffin, device = "pdf", width = 3, height = 8, units = "in")
