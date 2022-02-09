library(tidyverse)
library(ggfortify)
library(ggrepel)
library(umap)
library(lemon)
library(grid)

### Set working variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project/griffin/TFBS"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Uveal Melanoma/Manuscript/Figures/Figure 5"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project/TGL48_UVM_sample_list.txt"

### Get the list of files
filenames <- list.files(path, pattern = "*features*", recursive = TRUE, full.names = TRUE)
  
### Read in data and rbind into one dataframe
samples <- read.delim(samples)
datalist <- lapply(filenames, function(x){read.delim(file = x, header = TRUE)})
griffin <- do.call(rbind, datalist)
griffin <- as.data.frame(t(griffin))

### Make sure sample sheet and data have the same samples
samples <- samples[samples$sWGS %in% row.names(griffin), ]
colnames(griffin) <- griffin[1,]
griffin <- griffin[row.names(griffin) %in% samples$sWGS, ]
griffin[] <- mutate_all(griffin, function(x) as.numeric(as.character(x)))

### Format sample sheet
samples$Relapse <- factor(samples$Relapse, levels = c("Yes", "No", ""),
                          labels = c("Yes", "No", "LtFU"))
samples$Timepoint <- factor(samples$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))

### Format midpoint data
griffin_mid <- griffin[ , grep("midpoint", colnames(griffin))]
griffin_mid <- as.data.frame(scale(griffin_mid))

### Format amp data
griffin_amp <- griffin[ , grep("amp", colnames(griffin))]
griffin_amp <- as.data.frame(scale(griffin_amp))

### Format coverage data
griffin_cov <- griffin[ , grep("coverage", colnames(griffin))]
griffin_cov <- as.data.frame(scale(griffin_cov))

### Merge features together and order dataframe
griffin_scaled <- merge(griffin_mid, griffin_amp, by = "row.names")
griffin_scaled <- merge(griffin_scaled, griffin_cov, by.x = "Row.names", by.y = "row.names")
row.names(griffin_scaled) <- griffin_scaled$Row.names
griffin_scaled <- griffin_scaled[, -1]
griffin_scaled <- scale(griffin_scaled)
griffin_scaled <- griffin_scaled[samples$sWGS, ]

### Run UMAP
config <- umap.defaults
config$n_neighbors <- 4
config$negative_sample_rate <- 14

umap <- umap(griffin_scaled, config = config)
umap_df <- umap$layout %>% as.data.frame() %>%
  dplyr::rename(UMAP1="V1",
                UMAP2="V2") %>%
  dplyr::mutate(ID=row_number())
umap_df <- merge(umap_df, samples, by.x = "row.names", by.y = "sWGS")

### Check what it looks like
griffin <- ggplot(umap_df) + 
  geom_point(aes(UMAP1, UMAP2, color = Relapse, alpha = Timepoint)) +
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
  scale_color_manual(values = c("#fa4e0c", "#3f9e80", "#516eb0")) + 
  scale_alpha_manual(values = c(0.3, 0.5, 0.7, 0.9, 1))
griffin

write.table(umap_df, file.path(outdir, "Griffin_UMAP_table.txt"), sep = "\t")
