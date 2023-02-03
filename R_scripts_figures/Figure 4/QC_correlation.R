rm(list = ls())
cat("\014")

# Calcualte the correlation between QC metrics and methylation scores
library(ggplot2)
library("ggpubr")

setwd('')

qc_stats <- read.delim('Input_data/QC_metrics.tsv')

target <- c('02_Tumor','03_Tumor','05_Tumor','07_Tumor','08_Tumor','09_Tumor','10_Tumor',
            "03_BL","03_2W","03_3M","03_6M",
            "04_2W","04_3M","04_6M","04_12M",
            "05_BL","05_2W","05_3M","05_6M",
            "06_BL","06_2W","06_3M","06_6M","06_12M",
            "09_BL","09_2W","09_3M","09_6M",
            "01_BL","01_2W","01_3M","01_6M","01_12M",
            "02_BL","02_2W","02_3M","02_6M","02_12M",
            "11_BL","11_2W","11_3M","11_6M","11_12M",
            "07_BL","07_2W","07_3M","07_6M",
            "08_BL",
            "10_BL","10_2W","10_3M","10_6M","10_12M",
            "HBC-01-0015","HBC-01-0011","HBC-01-0014","HBC-01-0013","HBC-01-0012","HBC-01-0010",
            "HBC-01-0009","HBC-01-0017",
            "HBC-01-0001","HBC-01-0002","HBC-01-0003","HBC-01-0004","HBC-01-0005","HBC-01-0006")

qc_stats <- qc_stats[match(target, qc_stats$samplename),]

# cfmedip score
df <- read.delim('Input_data/methylation_score_cfmedip_0617.tsv')
score_HBC <- df[47:60,]
score_cfmedip <- df[1:46,]

score_medip <- read.delim('med2/tumor_score_0527.tsv')
score_medip <- score_medip[c("TGL48_0002_nn_P_PE_327_CM",
                             "TGL48_0003_nn_P_PE_328_CM",
                             "TGL48_0005_nn_P_PE_327_CM",
                             "TGL48_0007_nn_P_PE_312_CM",
                             "TGL48_0008_nn_P_PE_327_CM",
                             "TGL48_0009_nn_P_PE_325_CM",
                             "TGL48_0010_nn_P_PE_327_CM"),]


cfmedip_data <- data.frame(qc_stats[8:53,"Successfully.aligned.reads"])
cfmedip_data$scores <- score_cfmedip
colnames(cfmedip_data) <- c('reads','scores')

medip_data <- data.frame(qc_stats[1:7,"Successfully.aligned.reads"])
medip_data$scores <- score_medip
colnames(medip_data) <- c('reads','scores')

HBC_data <- data.frame(qc_stats[54:67,"Successfully.aligned.reads"])
HBC_data$scores <- score_HBC
colnames(HBC_data) <- c('reads','scores')

data <- rbind(medip_data, cfmedip_data)
data <- rbind(data, HBC_data)

data$Group <- c(rep('Tumor',7), rep('UMB',46), rep('HBC',14))

col.vec <- c(rep('#d73027',7), rep('#D55E00',46), rep('#91cf60',14))

pdf(file = "QC/QC_uveal_cor.pdf", width = 12, height = 12)
ggscatter(data, x = "reads", y = "scores",
          fill = "Group",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          palette = c("#91cf60","red","#D55E00"),
          xlab = "Successfully aligned reads (millions)", ylab = "Methylation score")
dev.off()
ggsave("QC_uveal_cor.pdf", width = 12, height = 12)
