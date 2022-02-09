library(tidyverse)
library(dplyr)

### Set working variables
path <- "/Users/derekwong/Desktop/H4H/projects/uveal_melanoma/insert_size/output"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project/sample_list.txt"
mouliere <- "/Users/derekwong/OneDrive - UHN/Post-Doc/External_data/Mouliere_fragment/Mouliere_fragment.txt"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project/insert_size"
project <- "TGL48_UVM"

swgs <- file.path(path, "swgs")

### Get file lists and names
insert_swgs <- list.files(path = swgs, pattern = "picard.txt", full.names = TRUE)

names_swgs <- list.files(path = swgs, pattern = "picard.txt", full.names = FALSE)
names_swgs <- gsub("(_WG).*", "\\1", names_swgs)

### Make outdir
dir.create(outdir, showWarnings = FALSE)

### Read in data
samples <- read.delim(samples)
mouliere <- read.delim(mouliere)

datalist <- lapply(insert_swgs, function(x){read.delim(file = x, skip = 13, colClasses = c(rep("integer", 2), rep("NULL", 2)), 
                                                       header = TRUE, col.names = c("length", "freq", "X", "Y"))})
insert_swgs <- Reduce(function(x,y) {merge(x,y, by = "length", all = TRUE, sort = TRUE)}, datalist)

### Subset only cfDNA samples
samples <- samples[samples$Type == "plasma", ]

### Format fragment data and save
row.names(insert_swgs) <- insert_swgs$length
insert_swgs <- insert_swgs[,-1]
colnames(insert_swgs) <- names_swgs
lengths <- c(10:600)
insert_swgs <- insert_swgs[row.names(insert_swgs) %in% lengths, colnames(insert_swgs) %in% samples$sWGS]
insert_swgs$length <- row.names(insert_swgs)
insert_swgs <- insert_swgs %>% dplyr::select(length, everything())
write.table(insert_swgs, file.path(outdir, paste0(project, "_fragment_swgs.txt")), row.names = FALSE, sep = "\t")

### Calculate fragment frequencies
freq <- insert_swgs[, -1]
sums <- colSums2(as.matrix(freq[row.names(freq) %in% c(10:250), ]))
freq <- as.data.frame(t(t(freq)/sums*100))
freq$length <- row.names(freq)
freq <- freq %>% dplyr::select(length, everything())
write.table(freq, file.path(outdir, paste0(project, "_fragment_swgs_freq.txt")), row.names = FALSE, sep = "\t")

### Calculate fragment proportions
short <- freq[ , -1]
short <- short[rownames(short) %in% c(20:150), ]
prop <- colSums2(as.matrix(short))/100
prop <- as.data.frame(prop)

### Format fragment proportions and save
prop$sample <- colnames(short)
prop$type <- "uveal_melanoma"
prop$diagnosis <- "uveal_melanoma"
colnames(prop) <- c("P.20_150.", "sample", "cancer", "type")
prop <- prop[ , c("sample", "cancer", "type", "P.20_150.")]
prop <- rbind(prop, mouliere)
write.table(prop, file.path(outdir, paste0(project, "_proportions.txt")), row.names = FALSE, sep = "\t")
