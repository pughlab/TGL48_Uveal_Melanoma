rm(list = ls())
cat("\014")

# Draw QC plot
library(ggplot2)

setwd('/Users/pluo/Project/UVM/QC')

files <- list.files('Input_data/QC/', recursive=F, full.names=F)

sample.medip <- c()
sample.cfmedip <- c()

total.reads.cfmedip <- c()
aligned.reads.cfmedip <- c()
methyl_F19K16.cfmedip <- c()
frac.duplicates.cfmedip <- c()
fracwoCpG.cfmedip <- c()
enrichrelH.cfmedip <- c()
enrichGoGe.cfmedip <- c()
deduplicates.cfmedip <- c()

aligned.reads.medip <- c()
methyl_F19K16.medip <- c()
frac.duplicates.medip <- c()
fracwoCpG.medip <- c()
enrichrelH.medip <- c()
enrichGoGe.medip <- c()
deduplicates.medip <- c()


for (f in files) {
  tmp <- strsplit(f, '_')[[1]]
  if (tmp[3] == 'nn') {
    df <- read.delim(paste0('Input_data/QC/',f))
    sample.medip <- c(sample.medip, colnames(df)[3])
    aligned.reads.medip <- c(aligned.reads.medip, df[df$metrics == 'ALIGNED_READS',3])
    methyl_F19K16.medip <- c(methyl_F19K16.medip, df[df$metrics == 'Pctg_BACs_is_methyl_F19K16',3])
    frac.duplicates.medip <- c(frac.duplicates.medip, as.numeric(df[df$metrics == 'READ_PAIR_DUPLICATES',3])/as.numeric(df[df$metrics == 'READ_PAIRS_EXAMINED',3]))
    fracwoCpG.medip <- c(fracwoCpG.medip, df[df$metrics == 'CpG_cov.fracReadsWoCpG',3])
    enrichrelH.medip <- c(enrichrelH.medip, df[df$metrics == 'enrich.enrichment.score.relH',3])
    enrichGoGe.medip <- c(enrichGoGe.medip, df[df$metrics == 'enrich.enrichment.score.GoGe',3])
    deduplicates.medip <- c(deduplicates.medip, as.numeric(df[df$metrics == 'READ_PAIRS_EXAMINED',3])-as.numeric(df[df$metrics == 'READ_PAIR_DUPLICATES',3]))
  } else {
    df <- read.delim(paste0('Input_data/QC/',f))
    sample.cfmedip <- c(sample.cfmedip, colnames(df)[3])
    #total.reads.cfmedip <- c(total.reads.cfmedip, df[df$QC_type == 'picard.CollectQualityYieldMetrics' & df$metrics == 'TOTAL_READS',3])
    total.reads.cfmedip <- c(total.reads.cfmedip, df[df$metrics == 'Total_reads',3])
    aligned.reads.cfmedip <- c(aligned.reads.cfmedip, df[df$metrics == 'ALIGNED_READS',3])
    methyl_F19K16.cfmedip <- c(methyl_F19K16.cfmedip, df[df$metrics == 'Pctg_BACs_is_methyl_F19K16',3])
    frac.duplicates.cfmedip <- c(frac.duplicates.cfmedip, as.numeric(df[df$metrics == 'READ_PAIR_DUPLICATES',3])/as.numeric(df[df$metrics == 'READ_PAIRS_EXAMINED',3]))
    fracwoCpG.cfmedip <- c(fracwoCpG.cfmedip, df[df$metrics == 'CpG_cov.fracReadsWoCpG',3])
    enrichrelH.cfmedip <- c(enrichrelH.cfmedip, df[df$metrics == 'enrich.enrichment.score.relH',3])
    enrichGoGe.cfmedip <- c(enrichGoGe.cfmedip, df[df$metrics == 'enrich.enrichment.score.GoGe',3])
    deduplicates.cfmedip <- c(deduplicates.cfmedip, as.numeric(df[df$metrics == 'READ_PAIRS_EXAMINED',3])-as.numeric(df[df$metrics == 'READ_PAIR_DUPLICATES',3]))
  }
}

# MeDIP
alignstats.medip <- data.frame(c('02_Tumor','03_Tumor','05_Tumor','07_Tumor','08_Tumor','09_Tumor','10_Tumor'))
colnames(alignstats.medip) <- 'samplename'
alignstats.medip$Successfully.aligned.reads <- as.numeric(as.character(aligned.reads.medip))/1000000
alignstats.medip$Pctg_BACs_is_methyl_F19K16 <- as.numeric(as.character(methyl_F19K16.medip))/100
alignstats.medip$duplicate.rates <- as.numeric(as.character(frac.duplicates.medip))
alignstats.medip$CpG_cov.fracReadsWoCpG <- as.numeric(as.character(fracwoCpG.medip))
alignstats.medip$enrich.enrichment.score.relH <- as.numeric(as.character(enrichrelH.medip))
alignstats.medip$enrich.enrichment.score.GoGe <- as.numeric(as.character(enrichGoGe.medip))
alignstats.medip$READ_PAIR_DEDUPLICATES <- as.numeric(as.character(deduplicates.medip))/1000000

# cfMeDIP
sample.name <- read.delim('Input_data/sample_name.tsv')$x

alignstats.cfmedip <- data.frame(sample.name)
colnames(alignstats.cfmedip) <- 'samplename'
#alignstats.cfmedip$Total.sequencing.reads <- as.numeric(as.character(total.reads.cfmedip))/1000000
alignstats.cfmedip$Successfully.aligned.reads <- as.numeric(as.character(aligned.reads.cfmedip))/1000000
alignstats.cfmedip$Pctg_BACs_is_methyl_F19K16 <- as.numeric(as.character(methyl_F19K16.cfmedip))/100
alignstats.cfmedip$duplicate.rates <- as.numeric(as.character(frac.duplicates.cfmedip))
alignstats.cfmedip$CpG_cov.fracReadsWoCpG <- as.numeric(as.character(fracwoCpG.cfmedip))
alignstats.cfmedip$enrich.enrichment.score.relH <- as.numeric(as.character(enrichrelH.cfmedip))
alignstats.cfmedip$enrich.enrichment.score.GoGe <- as.numeric(as.character(enrichGoGe.cfmedip))
alignstats.cfmedip$READ_PAIR_DEDUPLICATES <- as.numeric(as.character(deduplicates.cfmedip))/1000000

#names(alignstats.cfmedip)[names(alignstats.cfmedip) == 'sample.name'] <- 'samplename'
#write.table(alignstats.cfmedip$samplename,'sample_name.tsv',quote=F,row.names = F)

# sort based on total sequencing reads
# alignstats.cfmedip$samplename <- factor(alignstats.cfmedip$samplename,
#                                 levels = alignstats.cfmedip$samplename[order(alignstats.cfmedip$Total.sequencing.reads)])

# -------------- HBC -------------- #
files <- list.files('Input_data/QC_HBC/', recursive=F, full.names=F)
sample.HBC <- c()

aligned.reads.HBC <- c()
methyl_F19K16.HBC <- c()
frac.duplicates.HBC <- c()
fracwoCpG.HBC<- c()
enrichrelH.HBC <- c()
enrichGoGe.HBC <- c()
deduplicates.HBC <- c()

for (f in files) {
  df <- read.delim(paste0('Input_data/QC_HBC/',f))
  sample.HBC <- c(sample.HBC, colnames(df)[3])
  aligned.reads.HBC <- c(aligned.reads.HBC, df[df$metrics == 'ALIGNED_READS',3])
  methyl_F19K16.HBC <- c(methyl_F19K16.HBC, df[df$metrics == 'Pctg_BACs_is_methyl_F19K16',3])
  frac.duplicates.HBC <- c(frac.duplicates.HBC, as.numeric(df[df$metrics == 'READ_PAIR_DUPLICATES',3])/as.numeric(df[df$metrics == 'READ_PAIRS_EXAMINED',3]))
  fracwoCpG.HBC <- c(fracwoCpG.HBC, df[df$metrics == 'CpG_cov.fracReadsWoCpG',3])
  enrichrelH.HBC <- c(enrichrelH.HBC, df[df$metrics == 'enrich.enrichment.score.relH',3])
  enrichGoGe.HBC <- c(enrichGoGe.HBC, df[df$metrics == 'enrich.enrichment.score.GoGe',3])
  deduplicates.HBC <- c(deduplicates.HBC, as.numeric(df[df$metrics == 'READ_PAIRS_EXAMINED',3])-as.numeric(df[df$metrics == 'READ_PAIR_DUPLICATES',3]))
}

alignstats.HBC <- data.frame(c('HBC-01-0001','HBC-01-0002','HBC-01-0003','HBC-01-0004','HBC-01-0005','HBC-01-0006',
                               'HBC-01-0015','HBC-01-0011','HBC-01-0014','HBC-01-0013','HBC-01-0012','HBC-01-0010',
                               'HBC-01-0009','HBC-01-0017'))
colnames(alignstats.HBC) <- 'samplename'
alignstats.HBC$Successfully.aligned.reads <- as.numeric(as.character(aligned.reads.HBC))/1000000
alignstats.HBC$Pctg_BACs_is_methyl_F19K16 <- as.numeric(as.character(methyl_F19K16.HBC))/100
alignstats.HBC$duplicate.rates <- as.numeric(as.character(frac.duplicates.HBC))
alignstats.HBC$CpG_cov.fracReadsWoCpG <- as.numeric(as.character(fracwoCpG.HBC))
alignstats.HBC$enrich.enrichment.score.relH <- as.numeric(as.character(enrichrelH.HBC))
alignstats.HBC$enrich.enrichment.score.GoGe <- as.numeric(as.character(enrichGoGe.HBC))
alignstats.HBC$READ_PAIR_DEDUPLICATES <- as.numeric(as.character(deduplicates.HBC))/1000000

# sort based on sample
target <- c('02_Tumor','03_Tumor','05_Tumor','07_Tumor','08_Tumor','09_Tumor','10_Tumor',
            "01_BL","01_2W","01_3M","01_6M","01_12M",
            "02_BL","02_2W","02_3M","02_6M","02_12M",
            "03_BL","03_2W","03_3M","03_6M",
            "04_2W","04_3M","04_6M","04_12M",
            "05_BL","05_2W","05_3M","05_6M",
            "06_BL","06_2W","06_3M","06_6M","06_12M",
            "07_BL","07_2W","07_3M","07_6M",
            "08_BL",
            "09_BL","09_2W","09_3M","09_6M",
            "10_BL","10_2W","10_3M","10_6M","10_12M",
            "11_BL","11_2W","11_3M","11_6M","11_12M",
            'HBC-01-0001','HBC-01-0002','HBC-01-0003','HBC-01-0004','HBC-01-0005','HBC-01-0006',
            'HBC-01-0009','HBC-01-0010','HBC-01-0011','HBC-01-0012','HBC-01-0013','HBC-01-0014',
            'HBC-01-0015','HBC-01-0017')

#alignstats.cfmedip <- alignstats.cfmedip[match(target, alignstats.cfmedip$samplename),]
alignstats <- rbind(alignstats.medip, alignstats.cfmedip)
alignstats <- rbind(alignstats, alignstats.HBC)

write.table(alignstats, 'QC_metrics.tsv', row.names=F, quote = F, sep='\t')

alignstats$samplename <- factor(alignstats$samplename, levels=target)

# -------------- start to draw -------------- #

align_plots1 <- function (...) {
  pl <- list(...)
  stopifnot(do.call(all, lapply(pl, inherits, "gg")))
  gl <- lapply(pl, ggplotGrob)
  bind2 <- function(x, y) gtable:::rbind_gtable(x, y, "first")
  combined <- Reduce(bind2, gl[-1], gl[[1]])
  wl <- lapply(gl, "[[", "widths")
  combined$widths <- do.call(grid::unit.pmax, wl)
  grid::grid.newpage()
  grid::grid.draw(combined)
}

mytheme <- theme(axis.title.y = element_text(size = 25),
                 axis.title.x = element_blank(),
                 axis.line = element_line(color = "black"),
                 axis.text = element_text(size = 22),
                 axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "white", colour = "white"),
        plot.margin = unit(c(.5,0,0,0),"cm"))


# plot total sequencing reads plot
# myplot_totalseq <- ggplot(data = alignstats.cfmedip, aes(x = samplename)) +
#   geom_point(aes(y = Total.sequencing.reads), color = "black",size = 5) +
#   mytheme + theme(axis.text.x = element_blank()) +
#   scale_y_continuous(limits = c(0, max(alignstats.cfmedip$Total.sequencing.reads)),
#                      sec.axis = dup_axis(name = "Aligned reads"))

# barplot percentage aligned reads
# myplot_percaligned <- ggplot(aes(x = samplename,
#                                  y = Successfully.aligned.reads/Total.sequencing.reads),
#                              data = alignstats.cfmedip) +
#   geom_bar(stat = "identity") + mytheme + theme(axis.text.x = element_blank())

# barplot duplicated reads
# myplot_duplicates <- ggplot(aes(x = samplename,
#                                  y = duplicate.rates),
#                              data = alignstats) +
#   scale_y_continuous(limits = c(0, 1)) +
#   geom_bar(stat = "identity") + mytheme + theme(axis.text.x = element_blank())

# plot deduplicate sequencing reads plot
col.vec <- c(rep('#d73027',7), rep('#D55E00',46), rep('#91cf60',14))
myplot_deduplicates <- ggplot(aes(x = samplename, y = READ_PAIR_DEDUPLICATES),
                             data = alignstats) +
  geom_bar(stat = "identity", fill=col.vec) +
  mytheme + theme(axis.text.x = element_blank())

# barplot relH CpG enrichment score
myplot_relHScore <- ggplot(aes(x = samplename,
                               y = enrich.enrichment.score.relH),
                           data = alignstats) +
  #scale_y_continuous(limits = c(0, 1)) +
  geom_bar(stat = "identity", fill=col.vec) + mytheme + theme(axis.text.x = element_blank())
# barplot GoGe CpG enrichment score
myplot_GoGeScore <- ggplot(aes(x = samplename,
                               y = enrich.enrichment.score.GoGe),
                           data = alignstats) +
  #scale_y_continuous(limits = c(0, 1)) +
  geom_bar(stat = "identity", fill=col.vec) + mytheme + theme(axis.text.x = element_blank())

# barplot frac of read w/o CpG
myplot_fracWoCpG <- ggplot(aes(x = samplename,
                                y = CpG_cov.fracReadsWoCpG),
                            data = alignstats) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_bar(stat = "identity", fill=col.vec) + mytheme + theme(axis.text.x = element_blank())

# Arabidopsis
#Pctg_BACs_is_methyl_F19K16
#Pctg_BACs_is_unmethyl_F24B22
myplot_Arabidopsis <- ggplot(aes(x = samplename,
                                 y = Pctg_BACs_is_methyl_F19K16),
                             data = alignstats) +
  geom_bar(stat = "identity", fill=col.vec) + mytheme

pdf(file = "Suppl_Figure_4B.pdf", width = 12, height = 25, useDingbats = FALSE, onefile = FALSE)
align_plots1(
  # myplot_totalseq + #Add aligned reads as secondary y axis
  # geom_point(aes(y = Successfully.aligned.reads), color = "red",size = 5) +
  # theme(axis.title.y.right = element_text(color = "red"))+
  # ylab("Total sequencing reads (millions)"),
  # myplot_percaligned + ylab("Percentage\naligned reads"),
             # myplot_duplicates + ylab("Percentage\nduplicate reads"),
             myplot_deduplicates + ylab("Deduplicate\nread pairs (millions)"),
             myplot_relHScore + ylab("relH\nenrichment score"),
             myplot_GoGeScore + ylab("GoGe\nenrichment score"),
             myplot_fracWoCpG + ylab("Percentage\nreads without CpG"),
             myplot_Arabidopsis + ylab("Percentage\nspike-in reads methylated")
             )
dev.off()
#ggsave('QC_plot.pdf', re, device = "pdf", width = 15, height = 20, units = "in")

#test <- alignstats[alignstats$Pctg_BACs_is_methyl_F19K16 < 0.95,]
# min(alignstats$Pctg_BACs_is_methyl_F19K16)
# max(alignstats$Pctg_BACs_is_methyl_F19K16)
# min(alignstats$enrich.enrichment.score.relH)
# min(alignstats$enrich.enrichment.score.GoGe)
#
# max(alignstats$enrich.enrichment.score.relH)
# max(alignstats$enrich.enrichment.score.GoGe)
#
# #test <- alignstats.medip[alignstats.medip$READ_PAIR_DEDUPLICATES < 10000000]
# mean(alignstats$READ_PAIR_DEDUPLICATES)
# sd(alignstats$READ_PAIR_DEDUPLICATES)
