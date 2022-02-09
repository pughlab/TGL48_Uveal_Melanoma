library(tidyverse)
library(plyr)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/uveal_melanoma/fragmentomics/output"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project/sample_list.txt"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project/fragmentomics"
project <- "TGL48_UVM"

filenames_100kb <- list.files(path = path, pattern = "100kb_bins.txt", recursive = TRUE, full.names = TRUE)
filenames_5mb <- list.files(path = path, pattern = "5Mb_bins.txt", recursive = TRUE, full.names = TRUE)
filenames_dist <- list.files(path = path, pattern = "5Mb_dist.txt", recursive = TRUE, full.names = TRUE)
filenames_sum <- list.files(path = path, pattern = "summary.txt", recursive = TRUE, full.names = TRUE)
names <- list.files(path = path, pattern = "5Mb_bins.txt", recursive = TRUE, full.names = FALSE)
names <- sub("(_WG).*", '\\1', names)

### Make outdir
dir.create(outdir, showWarnings = FALSE)

### Read in data
samples <- read.delim(samples)

rows_100kb <- filenames_100kb[1]
rows_100kb <- read.delim(rows_100kb, header = TRUE)[ , 1:3]

rows_5mb <- filenames_5mb[1]
rows_5mb <- read.delim(rows_5mb, header = TRUE)[ , 2:5]

datalist <- lapply(filenames_100kb, function(x){read.delim(file = x, header = TRUE)[, 17]})
ratio_100kb <- do.call(cbind, datalist)
ratio_100kb <- as.data.frame(ratio_100kb)

datalist <- lapply(filenames_5mb, function(x){read.delim(file = x, header = TRUE)[, 16]})
ratio_5mb <- do.call(cbind, datalist)
ratio_5mb <- as.data.frame(ratio_5mb)

datalist <- lapply(filenames_dist, function(x){read.delim(file = x, header = TRUE)[, 9]})
dist <- do.call(cbind, datalist)
dist <- as.data.frame(dist)

summaries <- lapply(filenames_sum, function(x){read.delim(file = x, header = TRUE)})
summaries <- lapply(summaries, function(x){as.vector(t(x))})
summaries <- do.call(rbind, summaries)
summaries <- as.data.frame(summaries)

### Format fragmentation ratios and write data
samples <- samples[samples$Type == "plasma", ]
colnames(ratio_100kb) <- names
ratio_100kb <- ratio_100kb[ , colnames(ratio_100kb) %in% samples$sWGS]
ratio_100kb <- cbind(rows_100kb, ratio_100kb)
write.table(ratio_100kb, file.path(outdir, paste0(project, "_fragment_ratio_100kb.txt")), row.names = FALSE, sep = "\t")

colnames(ratio_5mb) <- names
ratio_5mb <- ratio_5mb[ , colnames(ratio_5mb) %in% samples$sWGS]
ratio_5mb <- cbind(rows_5mb, ratio_5mb)
write.table(ratio_5mb, file.path(outdir, paste0(project, "_fragment_ratio_5Mb.txt")), row.names = FALSE, sep = "\t")

### Format fragmentation distances and write data
colnames(dist) <- names
dist <- dist[ , colnames(dist) %in% samples$sWGS]
dist <- cbind(rows_5mb, dist)
write.table(dist, file.path(outdir, paste0(project, "_fragment_dist.txt")), row.names = FALSE, sep = "\t")

### Format fragmentomics summaries and write data
summaries <- summaries[ , c(2,7,8,13,14,19,21,23,25)]
summaries$V2 <- sub("(_WG).*", '\\1', summaries$V2)
summaries <- summaries[summaries$V2 %in% samples$sWGS, colSums(is.na(summaries))<nrow(summaries)]
colnames(summaries) <- c("sample", "ratio_cor", "ratio_sd", "coverage_cor", "coverage_dist", "nfrags", "mode", "mean", "median")
sum_names <- summaries$sample
summaries <- lapply(summaries, as.numeric)
summaries$sample <- sum_names
write.table(summaries, file.path(outdir, paste0(project, "_summary.txt")), row.names = FALSE, sep = "\t")
