library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(ggpubr)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project"
healthy_dir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Healthy_control_cohorts/combined_cohort"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Uveal Melanoma/Manuscript/Figures/Figure 2"

data_proportion <- file.path(path, "insert_size", "TGL48_UVM_proportions.txt")
data_samples <- file.path(path, "TGL48_UVM_sample_list.txt")

healthy_proportions <- file.path(healthy_dir, "insert_size", "HBC_fragment_proportions.txt")

### Read in data
data_proportion <- read.delim(data_proportion)
data_samples <- read.delim(data_samples)
healthy_proportions <- read.delim(healthy_proportions)

### Merge cohort fragment proportions with healthy control cohort
data_uveal <- data_proportion[data_proportion$type == "uveal_melanoma",]

### Format sample sheet
data_samples$Relapse <- factor(data_samples$Relapse, levels = c("No", "Yes", ""),
                               labels = c("No", "Yes", "Lost to Follow-up"))
data_samples$detection <- ifelse(data_samples$TF >= 0.03 | data_samples$TF_short >= 0.03, "Yes", "No")

### Format fragment proportions dataframe and set colors
colnames(data_proportion) <- c("sample", "cancer", "type", "p_20_150")
data_proportion$cancer <- factor(data_proportion$cancer,
                                 levels = c("healthy",  "renal", "glioblastoma", "bladder", "pancreatic", "melanoma", "breast", 
                                            "ovarian", "uveal_melanoma", "lung", "colorectal", "cholangiocarcinoma", "hepatocellular"),
                                 labels = c("Healthy", "KIRP", "GBM", "BLCA", "PAAD", "SKCM", "BRCA", "OVCA", "UVM", 
                                            "LUAD", "COAD", "CHOL", "LIHC"))
data_proportion <- data_proportion[complete.cases(data_proportion), ]
colors <- c(healthy = "lightgrey", low = "#FFFFB3", uveal_melanoma = "#fb9a99", high = "#CAB2D6")

### Calculate p-values for fragment proportions
data_stats <- data_proportion %>%
  group_by(cancer)%>% 
  dplyr::summarise(Median=median(p_20_150, na.rm = TRUE),
            Mean=mean(p_20_150, na.rm = TRUE),
            SD=sd(p_20_150, na.rm = TRUE),
            N=n())

c <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="KIRP"])$p.value
d <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="GBM"])$p.value
e <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="BLCA"])$p.value
f <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="PAAD"])$p.value
g <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="SKCM"])$p.value
h <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="BRCA"])$p.value
i <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="OVCA"])$p.value
j <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="UVM"])$p.value
k <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="LUAD"])$p.value
l <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="COAD"])$p.value
m <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="CHOL"])$p.value
n <- t.test(data_proportion$p_20_150[data_proportion$cancer=="Healthy"],data_proportion$p_20_150[data_proportion$cancer=="LIHC"])$p.value
t_test <- c(1, c, d, e, f, g, h, i ,j ,k ,l, m, n)
data_stats$pvalue <- t_test
data_stats$annot <- ifelse(data_stats$pvalue < 0.05 & data_stats$pvalue > 0.01, "*",
                           ifelse(data_stats$pvalue < 0.01 & data_stats$pvalue > 0.001, "**",
                                  ifelse(data_stats$pvalue < 0.001, "***", "")))
rm(c,d,e,f,g,h,i,j,k,l,m,t_test)

### Graph Figure 2A - proportion of small fragments across cohorts and cancers
FigA <- ggplot(data_proportion, aes(cancer, p_20_150) ) +
  geom_boxplot(aes(fill = type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 2, alpha = 0.5, width = 0.15) +
  geom_text(data = data_stats, aes(x = cancer, y = 0, label = N), size = 4) +
  geom_text(data = data_stats, aes(x = cancer, y = 0.55, label = annot), size = 5) +
  scale_fill_manual(labels = c("Healthy", "Low ctDNA", "Uveal Melanoma", "High ctDNA"), values = colors) +
  labs(fill = "") +
  xlab("") + 
  ylab("Proportion under 150bp") +
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(face = c(replicate(8, "plain"), "bold", replicate(4, "plain")),
                                   angle = 0, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 13),
        legend.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  scale_y_continuous(limits=c(-0.05, 0.6), expand = c(0,0)) +
  guides(fill=guide_legend(nrow = 2,byrow = TRUE)) +
  coord_flip()
FigA
#ggsave(file.path(outdir, "Figure 2 - Fragment proportion.pdf"), FigA, device = "pdf", width = 5, height = 7, units = "in")

### Separate out uveal melanoma from fragment proportions and append clinical information
data_uveal <- rbind(data_uveal, healthy_proportions)
data_uveal <- merge(data_uveal, data_samples, by.x = "sample", by.y = "sWGS", all = TRUE)
data_uveal <- data_uveal[!(is.na(data_uveal$P.20_150.)), ]
data_uveal$Relapse <- ifelse(is.na(data_uveal$Relapse), "Healthy Control", data_uveal$Relapse)
data_uveal$Relapse <- factor(data_uveal$Relapse, levels = c("Healthy Control", 1, 2, 3),
                             labels = c("Healthy", "Remission", "Relapsed", "Lost to\nFollow-up"))

### Calculate stats and p-values for uveal cohort fragment proportions (Relapse)
data_stats_uveal <- data_uveal %>%
  group_by(Relapse)%>% 
  dplyr::summarise(Median=median(P.20_150., na.rm = TRUE),
                   Mean=mean(P.20_150., na.rm = TRUE),
                   SD=sd(P.20_150., na.rm = TRUE),
                   N=n())
a <- t.test(data_uveal$P.20_150.[data_uveal$Relapse=="Healthy"], data_uveal$P.20_150.[data_uveal$Relapse=="Remission"])$p.value
b <- t.test(data_uveal$P.20_150.[data_uveal$Relapse=="Healthy"], data_uveal$P.20_150.[data_uveal$Relapse=="Relapsed"])$p.value
c <- t.test(data_uveal$P.20_150.[data_uveal$Relapse=="Healthy"], data_uveal$P.20_150.[data_uveal$Relapse=="Lost to\nFollow-up"])$p.value
d <- t.test(data_uveal$P.20_150.[data_uveal$Relapse=="Remission"], data_uveal$P.20_150.[data_uveal$Relapse=="Relapsed"])$p.value
e <- t.test(data_uveal$P.20_150.[data_uveal$Relapse=="Remission"], data_uveal$P.20_150.[data_uveal$Relapse=="Lost to\nFollow-up"])$p.value
t_test <- c(1, a, b, c)
t_test2 <- c(1, 1, d, e)
data_stats_uveal$pvalue <- t_test
data_stats_uveal$annot <- ifelse(data_stats_uveal$pvalue < 0.05 & data_stats_uveal$pvalue > 0.01, "*",
                           ifelse(data_stats_uveal$pvalue < 0.01 & data_stats_uveal$pvalue > 0.001, "**",
                                  ifelse(data_stats_uveal$pvalue < 0.001, "***", "")))
rm(a,b,t_test)

## Graph Figure 2B - Fragment proportions within Uveal cohort (Relapse)
FigB <- ggplot(data_uveal, aes(Relapse, P.20_150.)) +
  geom_boxplot(aes(fill = Relapse), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 2, alpha = 0.5, width = 0.15) +
  geom_text(data = data_stats_uveal, aes(x = Relapse, y = 0, label = N), size = 4) +
  geom_text(data = data_stats_uveal, aes(x = Relapse, y = 0.28, label = annot), size = 5.5) +
  xlab("") + 
  ylab("Proportion under 150bp") +
  ggtitle("") + 
  scale_fill_manual(labels = c("Healthy Control", "No", "Yes", "Lost to\nFollow-up"), 
                    values = c("grey", "grey", "#fb9a99", "grey")) +
  labs(fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(face = c(replicate(2, "plain"), "bold", "plain"),
                                   angle = 0, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 1)) + 
  scale_y_continuous(limits=c(-0.03, 0.3), expand = c(0,0)) +
  guides(fill=guide_legend(nrow = 2,byrow = TRUE)) +
  coord_flip()
FigB
#ggsave(file.path(outdir, "Figure 2 - Fragment proportion cohort.pdf"), FigB, device = "pdf", width = 5, height = 2.5, units = "in")

Figure <- ggpubr::ggarrange(FigA, FigB, heights = c(8,4), nrow = 2, align = "v")
Figure
ggsave(file.path(outdir, "Figure 2 - Fragment proportion.pdf"), Figure, device = "pdf", width = 4, height = 9, units = "in")
