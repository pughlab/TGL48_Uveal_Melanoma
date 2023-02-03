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

eye <- file.path(path, "griffin/uveal/TGL48_UVM_griffin_corrected_Eye.txt")
samples <- file.path(path, "sample_list.txt")
healthy_eye <- file.path(healthy, "HBC_griffin_corrected_Eye.txt")

### Import data
eye <- read.delim(eye)
samples <- read.delim(samples)
healthy_eye <- read.delim(healthy_eye, check.names = FALSE)

### Format data
distance <- eye$distance
eye <- eye[ , -1]
healthy_eye <- healthy_eye[ , -1]

### Scale data
eye <- as.data.frame(scale(eye)*0.1 + 1)
healthy_eye <- as.data.frame(scale(healthy_eye)*0.1 + 1)

### Subset samples to compare (09 vs 11)
samples <- samples[samples$Patient %in% c("UMB-009", "UMB-011") &
                     !(samples$Timepoint %in% c("Lymphocytes", "Tumour")), ]
eye <- eye[, samples$sWGS]

eye$distance <- distance

### Make healthy median and sd
eye_median <- rowMedians(as.matrix(healthy_eye))
eye_sd <- rowSds(as.matrix(healthy_eye))
eye_median <- as.data.frame(cbind(eye_median, eye_sd))
colnames(eye_median) <- c("median", "sd")
eye_median$distance <- distance
eye_median$variable <- "Healthy"
eye_median$Patient <- "Healthy"
eye_median$type <- "Eye"
eye_median$Timepoint <- "Baseline"

median <- eye_median

### Melt data
eye <- reshape2::melt(eye, id.vars = c("distance"))
eye$type <- "Eye"

### Format melted data
data_melt <- eye
data_melt <- merge(data_melt, samples, by.x = c("variable"), by.y = c("sWGS"))
data_melt$Timepoint <- factor(data_melt$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"),
                              labels = c("Baseline", "2 Weeks", "3 Months", "6 Months", "12 Months"))
data_melt$Patient <- factor(data_melt$Patient, levels = c("UMB-009", "UMB-011"),
                            labels = c("09 - Recurred", "11 - Remission"))
data_melt <- data_melt[order(data_melt$Patient,
                             data_melt$type,
                             data_melt$Timepoint,
                             data_melt$distance), ]

### Make reference line
data_ref <- data_melt[data_melt$Timepoint == "Baseline", c(2, 4, 5)]
data_ref$value <- ifelse(data_ref$type == "Eye", 0.8, 0.8)
healthy_ref <- data_ref
healthy_ref$Patient <- "Healthy"
data_ref <- rbind(data_ref, healthy_ref)

### Graph Figure - cohort wide
griffin <- ggplot(data_melt, aes(distance, value)) +
  geom_line(aes(group = variable, color = Timepoint), alpha = 0.75, size = 0.55) +
  geom_line(data = median, aes(distance, median, group = variable, color = Timepoint), alpha = 0.75, size = 0.55) +
  geom_ribbon(data = median, aes(distance, y = median, ymin = (median - sd), ymax = (median + sd), group = variable), alpha = 0.5) + 
  geom_hline(data = data_ref, aes(yintercept = value), linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("Baseline" = "black", "2 Weeks" = "#EF3B2C", "3 Months" = "#CB181D", "6 Months" = "#A50F15", "12 Months" = "#67000D")) + 
  facet_grid(.~Patient, scales = "free") +
  facet_rep_grid(.~Patient) + 
  labs(fill = "") +
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
ggsave(file.path(outdir, paste0("Supplemental Figure 12.pdf")), griffin, device = "pdf", width = 7.5, height = 4, units = "in")
