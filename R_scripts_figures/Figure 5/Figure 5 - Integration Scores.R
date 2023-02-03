library(ggplot2)
library(reshape2)
library(plyr)
library(ggforce)
library(ggbreak)
library(cowplot)

path <- ""
outdir <- ""
samples <- "TGL48_UVM_sample_list.txt"

### Import data
samples <- read.delim(samples)
data_vaf <- read.delim(file.path(path, "variant_vafs.txt"))
data_ichor <- read.delim(file.path(path, "ichor_score.txt"))
data_insert <- read.delim(file.path(path, "insert_size.txt"))
data_ratio <- read.delim(file.path(path, "ratios.txt"))
data_PRC <- read.delim(file.path(path, "PRC1_scores.txt"))
data_griffin <- read.delim(file.path(path, "griffin_metrics.txt"))
data_relapse <- read.delim(file.path(outdir, "relapse.txt"))

### Format samples data
samples <- samples[!(samples$Timepoint %in% c("Healthy", "Tumour", "Lymphocytes")), ]
samples$Relapse <- factor(samples$Relapse, levels = c("Yes", "No", ""),
                          labels = c("Relapse", "Remission", "Lost to Follow-up"))
samples$Timepoint <- factor(samples$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"),
                            labels = c("Baseline", "2 Weeks", "3 Months", "6 Months", "12 Months"))
samples <- samples[, c("Patient", "Timepoint", "Relapse", "cfMeDIP", "sWGS", "Targeted.Panel")]

### Merge metrics with samples sheet
data <- merge(samples, data_vaf, by.x = "Targeted.Panel", by.y = "sample", all = TRUE)
data <- merge(data, data_ichor, by.x = "Targeted.Panel", by.y = "sample", all = TRUE)
data <- merge(data, data_insert, by.x = "sWGS", by.y = "sample", all = TRUE)
data <- merge(data, data_ratio, by.x = "sWGS", by.y = "sample", all = TRUE)
data <- merge(data, data_PRC, by.x = "sWGS", by.y = "sample", all = TRUE)
data <- merge(data, data_griffin, by.x = "sWGS", by.y = "sample", all = TRUE)
colnames(data) <- c("sWGS", "Targeted.Panel", "Patient", "Timepoint", "Relapse", "cfMeDIP", "vaf", "ichor", "ichor_limit", 
                    "insert", "insert_limit", "ratio", "ratio_limit", "PRC1", "PRC1_limit", "griffin", "griffin_limit")
limits <- data[, c("Patient", "Timepoint", "Relapse", "ichor_limit", "insert_limit", "ratio_limit", "PRC1_limit", "griffin_limit")]
data <- data[, c("Patient", "Timepoint", "Relapse", "vaf", "ichor", "insert", "ratio", "PRC1", "griffin")]

### Order samples by clinical outcome
data_relapse$sample_ID <- gsub("UMB-0", "", data_relapse$sample_ID)
data_relapse <- data_relapse[order(data_relapse$relapse,
                                   data_relapse$death,
                                   data_relapse$followup), ]
order <- data_relapse$sample_ID

data$Patient <- gsub("UMB-0", "", data$Patient)
data$Patient <- factor(data$Patient, levels = order)

limits$Patient <- gsub("UMB-0", "", limits$Patient)
limits$Patient <- factor(limits$Patient, levels = order)

### Melt data
data_melt <- melt(data, id = c("Patient", "Timepoint", "Relapse"))
data_melt_limits <- melt(limits, id = c("Patient", "Timepoint", "Relapse"))

### Format tables for comparison plotting
data_melt$variable <- factor(data_melt$variable, levels = c("vaf", "ichor", "insert", "ratio", "PRC1", "griffin"),
                             labels = c("Variant", "Copy\nNumber", "Insert Size", "Fragment\nGenome", "Fragment\nPRC1", "Fragment\nCoverage"))

data_melt_limits$variable <- factor(data_melt_limits$variable, levels = c("ichor_limit", "insert_limit", "ratio_limit", "PRC1_limit", "griffin_limit"),
                                    labels = c("Copy\nNumber", "Insert Size", "Fragment\nGenome", "Fragment\nPRC1", "Fragment\nCoverage"))

### Plot Metrics comparisons by Patient
time_plot <- ggplot(data_melt, aes(x = Timepoint, y = value, fill = Relapse, color = Relapse)) +
  geom_point(shape = 16, size = 2) +
  geom_line(aes(group = Patient)) +
  #geom_line(data = data_melt_limits, aes(group = variable, color = "red"), linetype = "solid", size = 0.75, alpha = 0.75) +
  facet_grid(variable~Patient, scales = "free", space = "free_x") +
  xlab("Timepoint") + 
  ylab("Score") +
  ggtitle("Individual Metrics") + 
  labs(fill = "Timepoint") +
  scale_fill_manual(values = c("Relapse" = "#e14b31", "Remission" = "#22a7f0", "Lost to Follow-up" = "#969696"), guide = "none") +
  scale_color_manual(values = c("Relapse" = "#e14b31", "Remission" = "#22a7f0", "Lost to Follow-up" = "#969696", "red" = "black"), guide = "none") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 15),
        axis.title.x = element_text(size = 0),
        strip.background = element_rect(color="white", fill="white"),
        strip.text = element_text(size = 12, color = "black", face = "bold", angle = 0)) +
  scale_y_continuous(limits=c(-0.15, 1.2), expand = c(0,0), breaks = c(0, 0.5, 1))
time_plot
ggsave(file.path(outdir, "Figure 7 - Analysis Metrics.pdf"), time_plot, device = "pdf", width = 10, height = 7.5, units = "in")

### Combine metrics (scaled)
data$ichor <- data$ichor*0.1
data$insert <- data$insert*0.225
data$ratio <- data$ratio*0.225
data$PRC1 <- data$PRC1*0.225
data$griffin <- data$griffin*0.225
data$combined <- rowSums2(as.matrix(data[ , c("ichor", "insert", "ratio", "PRC1", "griffin")]), na.rm = TRUE)
#limits$combined <- rowMeans2(as.matrix(limits[, 4:7]))

### Normalize to baseline
#data$combined <- data$combined/limits$combined
baseline <- data[, c("Patient", "Relapse")]
baseline$Timepoint <- "Baseline"
baseline$Timepoint <- ifelse(baseline$Patient == "04", "2 Weeks", "Baseline")
baseline <- merge(baseline, data, by = c("Patient", "Relapse", "Timepoint"))
data$combined <- data$combined/baseline$combined

### Plot Metric means per patient
patient_plot <- ggplot(data, aes(x = Timepoint, y = combined)) +
  geom_point(aes(color = Relapse), shape = 16, size = 2) +
  geom_line(aes(color = Relapse, group = Patient)) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", size = 0.75, alpha = 0.75) +
  facet_grid(.~Patient, scales = "free", space = "free") +
  xlab("Timepoint") + 
  ylab("Score") +
  ggtitle("Combined Metric") + 
  labs(color = "Patient Outcome") +
  scale_color_manual(values = c("Relapse" = "#e14b31", "Remission" = "#22a7f0", "Lost to Follow-up" = "#969696", "red" = "black"), 
                     limits = c("Relapse", "Remission", "Lost to Follow-up")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13),
        axis.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        strip.background = element_rect(color="white", fill="white"),
        strip.text = element_text(size = 12, color = "black", face = "bold", angle = 0)) + 
  scale_y_continuous(limits=c(0, 2), expand = c(0,0), breaks = c(0, 1.0, 2.0))
patient_plot
ggsave(file.path(outdir, "Figure 7 - Combined.pdf"), device = "pdf", width = 10, height = 3.5, units = "in")

Figure <- plot_grid(time_plot, patient_plot, labels = c('A','B'), label_size = 20,
                    align = 'v', axis = 'lr', ncol = 1, rel_heights = c(1.75, 1))
Figure
ggsave(file.path(outdir, "Figure 7.pdf"), Figure, device = "pdf", width = 9, height = 10, units = "in")
