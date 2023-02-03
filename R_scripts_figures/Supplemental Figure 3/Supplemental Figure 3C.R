library(tidyverse)
library(dplyr)
library(reshape2)
library(ggformula)
library(signal)
library(matrixStats)
library(lemon)
library(ggpubr)

### Set paths
path <- ""
outdir <- ""
healthy <- ""

site <- c("Cardiac", "Lymphoid", "Digestive", "Musculoskeletal", "Myeloid_erythroid", "Neural", "Pulmonary", "Vascular_endothelial")
samples <- file.path(path, "sample_list.txt")
healthy_samples <- file.path("HBC_sample_list.txt")

### Read in samples
samples <- read.delim(samples)
healthy_samples <- read.delim(healthy_samples, check.names = FALSE)

### prune healthy samples
#remove <- healthy_samples[healthy_samples$cohort == "Phallen", ]
#remove <- remove$sample

### Get list of files
files <- file.path(path, paste0("griffin/DHS/TGL48_UVM_griffin_corrected_", site, ".txt"))
healthy_files <- file.path(healthy, paste0("HBC_griffin_corrected_", site, ".txt"))

### Import data
organ <- lapply(files, function(x){read.delim(file = x, header = TRUE)})
healthy_organ <- lapply(healthy_files, function(x){read.delim(file = x, header = TRUE)})

### Get row and column names
distance <- organ[[1]][["distance"]]
names <- colnames(organ[[1]])
names <- names[!(names == "distance")]
names(organ) <- site

healthy_names <- colnames(healthy_organ[[1]])
healthy_names <- healthy_names[!(healthy_names == "distance")]
names(healthy_organ) <- site

### Format data
my_fun <- function(x) {
  names(x) <- site
  x <- x[ , -1]
  x <- as.data.frame(scale(x)*0.1 + 1)
  row.names(x) <- distance
  x <- x[row.names(x) %in% c("-30", "-15", "0", "15", "30"), ]
  x <- colMedians(as.matrix(x))
  x <- as.data.frame(x)
}

organs <- lapply(organ, my_fun)
organs <- do.call(cbind, organs)
colnames(organs) <- site
row.names(organs) <- names

healthy_organs <- lapply(healthy_organ, my_fun)
healthy_organs <- do.call(cbind, healthy_organs)
colnames(healthy_organs) <- site
row.names(healthy_organs) <- healthy_names
#healthy_organs <- healthy_organs[!(row.names(healthy_organs) %in% remove), ]

### Subset samples to only plasmas and order samples
samples <- samples[!(samples$Timepoint %in% c("Lymphocytes", "Tumour", "Healthy")), ]
organs <- organs[row.names(organs) %in% samples$sWGS, ]
organs <- organs[samples$sWGS, ]
organs <- merge(organs, samples[, c("sWGS", "Patient", "Timepoint")], by.x = "row.names", by.y = "sWGS")
organs <- organs[, -1]
data_melt <- melt(organs, id = c("Patient", "Timepoint"))
data_melt$Patient <- factor(data_melt$Patient, levels = c(paste0("UMB-0", c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11"))))
data_melt$Timepoint <- factor(data_melt$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))

### Make healthy median and sd
medians <- colMedians(as.matrix(healthy_organs))
sds <- colSds(as.matrix(healthy_organs))

medians <- do.call(rbind, replicate(nrow(organs), medians, simplify = FALSE))
colnames(medians) <- site
medians <- melt(medians)

sds <- do.call(rbind, replicate(nrow(organs), sds, simplify = FALSE))
colnames(sds) <- site
sds <- melt(sds)

### Bind data and healthy medians
data_melt <- cbind(data_melt, medians[, c("value")])
colnames(data_melt) <- c("Patient", "Timepoint", "variable", "value", "median")
data_melt$upper <- data_melt$median + 1.5*sds$value
data_melt$lower <- data_melt$median - 1.5*sds$value

### Format melted data
data_melt$Timepoint <- factor(data_melt$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"),
                              labels = c("Baseline", "2 Weeks", "3 Months", "6 Months", "12 Months"))
data_melt$variable <- factor(data_melt$variable, levels = site,
                             labels = c("Cardiac", "Lymphoid", "Digestive", "Musculoskeletal", "Myeloid", "Neural", "Pulmonary", "Endothelial"))
data_melt <- data_melt[order(data_melt$variable,
                             data_melt$Patient,
                             data_melt$Timepoint), ]

### Graph Figure
site_plot <- ggplot(data_melt, aes(Timepoint, value)) +
  geom_point(aes(color = Timepoint), alpha = 0.75, size = 2) +
  geom_hline(aes(yintercept = median), linetype = "dashed") +
  geom_hline(aes(yintercept = upper), linetype = "dashed", alpha = 0.5) +
  geom_hline(aes(yintercept = lower), linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("Baseline" = "black", "2 Weeks" = "#EF3B2C", "3 Months" = "#CB181D", "6 Months" = "#A50F15", "12 Months" = "#67000D")) + 
  facet_grid(variable ~ Patient, scales = "free", space = "free") +
  facet_rep_grid(variable ~ Patient) + 
  labs(color = "") +
  xlab("Timepoint") + 
  ylab("Mean Central Coverage (+/- 30bp)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
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
        strip.text = element_text(face = "bold", size = 8),
        panel.spacing = unit(-3, "mm"),
        panel.spacing.y = unit(-10, "mm")) +
  guides(color=guide_legend(nrow = 3,byrow = TRUE)) + 
  scale_y_continuous(limits=c(0.55, 1.35), expand = c(0,0))
site_plot

ggsave(file.path(outdir, paste0("Supplemental Figure 3C.pdf")), site_plot, device = "pdf", width = 8, height = 10, units = "in")
