library(tidyverse)
library(plyr)

### Set working variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project"
project <- "TGL48_UVM"

samples <- file.path(path, "sample_list.txt")
coverage <- file.path(path, "GATK_coverage", "TGL48_UVM_coverage_summary.txt")
ichorCNA <- file.path(path, "ichorCNA", "TGL48_UVM_ichorCNA_summary.txt")

### Read in files
samples <- read.delim(samples, check.names = FALSE)
coverage <- read.delim(coverage)
ichorCNA <- read.delim(ichorCNA)

### Merge data
samples <- merge(samples, coverage, by.x = c("Patient", "Timepoint", "Targeted Panel", "sWGS"), by.y = c("patient_ID", "timepoint", "ts", "swgs"), all = TRUE)
samples <- merge(samples, ichorCNA, by.x = c("sWGS"), by.y = c("sample"), all = TRUE)
samples <- samples[ , c("Patient", "Timepoint", "Stage", "Treatment", "Sex", "Volume", "Age", "Relapse", "cfMeDIP", "sWGS", "Targeted Panel",
                        "coverage_swgs", "coverage_short", "coverage_all_unique", "coverage_ts", "TF", "TF_short")]
samples <- samples[order(samples$Patient, 
                         factor(samples$Timepoint, levels = c("Tumour", "Lymphocytes", "Baseline", "2 weeks", "3 months", "6 months", "12 months"))), ]

write.table(samples, file.path(path, paste0(project, "_sample_list.txt")), row.names = FALSE, sep = "\t")
