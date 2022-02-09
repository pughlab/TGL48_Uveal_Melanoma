library(tidyverse)
library(ggfortify)
library(ggforce)
library(ggrepel)
library(lemon)
library(grid)

### Set working variables
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Uveal Melanoma/Manuscript/Figures/Figure 5"

### Read in UMAP table
umap_df <- read.delim(file.path(outdir, "Griffin_UMAP_table.txt"))

### Format UMAP table
umap_df$Relapse <- factor(umap_df$Relapse, levels = c("Yes", "No", "LtFU"),
                          labels = c("Relapse", "Remission", "Lost to Follow-up"))
umap_df$Timepoint <- factor(umap_df$Timepoint, levels = c("Baseline", "2 weeks", "3 months", "6 months", "12 months"))

### Plot UMAP
griffin <- ggplot(umap_df) + 
  geom_point(aes(UMAP1, UMAP2, color = Relapse)) +
  geom_ellipse(aes(x0 = 2, y0 = -4.5, a = 3, b = 1.8, angle = 40), color = "#fa4e0c", linetype = "dashed") +
  geom_ellipse(aes(x0 = -2.4, y0 = 1.1, a = 4.75, b = 1.8, angle = 20), color = "#3f9e80", linetype = "dashed") +
  geom_ellipse(aes(x0 = 2, y0 = 4, a = 1, b = 2.5, angle = 0), color = "#fa4e0c", linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) + 
  scale_color_manual(values = c("#fa4e0c", "#3f9e80", "#516eb0")) + 
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(-7.5, 7.5))
griffin
#ggsave(file.path(outdir, paste0("Figure 5 - Griffin UMAP_cohort.pdf")), griffin, device = "pdf", width = 4.5, height = 3.5, units = "in")

### Format for per patient plotting
umap_df <- umap_df[order(factor(umap_df$Relapse, levels = c("Relapse", "Remission", "Lost to Follow-up"))), ]
patient_order <- umap_df$Patient
patient_order <- unique(patient_order)
umap_df$Patient <- factor(umap_df$Patient, levels = c(" ", patient_order))

### Plot per patient
patient <- ggplot(umap_df, aes(UMAP1, UMAP2)) + 
  geom_point(aes(color = Relapse)) +
  geom_path(size = 0.5, arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches"))) + 
  geom_ellipse(aes(x0 = 2, y0 = -4.5, a = 3, b = 1.8, angle = 40), color = "#fa4e0c", linetype = "dashed") +
  geom_ellipse(aes(x0 = -2.4, y0 = 1.1, a = 4.75, b = 1.8, angle = 20), color = "#3f9e80", linetype = "dashed") +
  geom_ellipse(aes(x0 = 2, y0 = 4, a = 1, b = 2.5, angle = 0), color = "#fa4e0c", linetype = "dashed") +
  labs(color="Patient Outcome") +
  facet_wrap(.~Patient, scales = "fixed", nrow = 2, drop = FALSE) +
  facet_rep_wrap(.~Patient, scales = "fixed", nrow = 2, drop = FALSE) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(-5, "mm"),
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 10, face = "bold"),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.position = c(0.05, 0.8),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) + 
  scale_color_manual(values=c("Relapse" = "#D55E00", "Remission" = "#359B73", "Lost to Follow-up" = "#2271B2"))
patient
#ggsave(file.path(outdir, paste0("Figure 5 - Griffin UMAP_patient.pdf")), griffin, device = "pdf", width = 9, height = 4, units = "in")

### Edit the grobs
g <- ggplotGrob(patient)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-1-1", "strip-t-1-1")
# remove grobs
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]
## move axis closer to panel
g$layout[g$layout$name == "axis-l-1-1", c("l", "r")] = c(8, 8)
g$layout[g$layout$name == "axis-b-1-1", c("t", "b")] = c(14, 14)
grid.newpage()
grid.draw(g)

### Combine into Figure
Figure <- ggpubr::ggarrange(griffin, g, widths = c(1.5, 4), nrow = 1)
Figure
ggsave(file.path(outdir, "Figure 5 - Griffin.pdf"), Figure, device = "pdf", width = 12, height = 3.5, units = "in")
