library(tidyverse)
library(dplyr)

# Set working variables
all_path <- "/Users/derekwong/Desktop/H4H/projects/uveal_melanoma/pipeline_output/all_unique/"
sscs_path <- "/Users/derekwong/Desktop/H4H/projects/uveal_melanoma/pipeline_output/sscs/"
dcs_path <- "/Users/derekwong/Desktop/H4H/projects/uveal_melanoma/pipeline_output/dcs/"
outdir <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Uveal_Melanoma_Project/pipeline_analysis"
project <- "TGL48_UVM"
date <- Sys.Date()

# Make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Get Files
germline_all <- list.files(path = file.path(all_path, "HaplotypeCaller/CPSR"), pattern = "*.tsv", full.names = TRUE)

mutect2_all <- list.files(path = file.path(all_path, "MuTect2"), pattern = "*.tsv", full.names = TRUE)
mutect2_sscs <- list.files(path = file.path(sscs_path, "MuTect2"), pattern = "*.tsv", full.names = TRUE)
mutect2_dcs <- list.files(path = file.path(dcs_path, "MuTect2"), pattern = "*.tsv", full.names = TRUE)

strelka_all <- list.files(path = file.path(all_path, "Strelka"), pattern = "*.tsv", full.names = TRUE)
strelka_sscs <- list.files(path = file.path(sscs_path, "Strelka"), pattern = "*.tsv", full.names = TRUE)

pindel_all <- list.files(path = file.path(all_path, "Pindel"), pattern = "*.tsv", full.names = TRUE)

# Import data
germline_all <- read_tsv(germline_all)

mutect2_all <- read_tsv(mutect2_all)
mutect2_sscs <- read_tsv(mutect2_sscs)
mutect2_dcs <- read_tsv(mutect2_dcs)

strelka_all <- read_tsv(strelka_all)
strelka_sscs <- read_tsv(strelka_sscs)

pindel_all <- read_tsv(pindel_all)

# Parse mutations (Mutect2)
mutect2_all <- mutect2_all[mutect2_all$Hugo_Symbol %in% c("BAP1", "GNA11", "GNAQ", "SF3B1", "EIF1AX"), ]
mutect2_all <- mutect2_all[, colnames(mutect2_all) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                                             "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                                             "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG")]
mutect2_all <- mutect2_all[!grepl("benign", mutect2_all$SIFT) &
                               !grepl("benign", mutect2_all$PolyPhen) &
                               !grepl("benign", mutect2_all$CLIN_SIG) &
                               !(mutect2_all$Variant_Classification %in% c("3'UTR", "5'UTR", "3'Flank", "Intron", "5'Flank", "RNA", "Silent", "Splice_Region")), ]

mutect2_sscs <- mutect2_sscs[mutect2_sscs$Hugo_Symbol %in% c("BAP1", "GNA11", "GNAQ", "SF3B1", "EIF1AX"), ]
mutect2_sscs <- mutect2_sscs[, colnames(mutect2_sscs) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                                                "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                                                "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG")]
mutect2_sscs <- mutect2_sscs[!grepl("benign", mutect2_sscs$SIFT) &
                               !grepl("benign", mutect2_sscs$PolyPhen) &
                               !grepl("benign", mutect2_sscs$CLIN_SIG) &
                               !(mutect2_sscs$Variant_Classification %in% c("3'UTR", "5'UTR", "3'Flank", "Intron", "5'Flank", "RNA", "Silent", "Splice_Region")), ]
mutect2_sscs <- anti_join(mutect2_sscs, mutect2_all, by = c("Start_Position", "Tumor_Sample_Barcode", "HGVSc"))

mutect2_dcs <- mutect2_dcs[mutect2_dcs$Hugo_Symbol %in% c("BAP1", "GNA11", "GNAQ", "SF3B1", "EIF1AX"), ]
mutect2_dcs <- mutect2_dcs[, colnames(mutect2_dcs) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                                             "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                                             "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG")]
mutect2_dcs <- mutect2_dcs[!grepl("benign", mutect2_dcs$SIFT) &
                               !grepl("benign", mutect2_dcs$PolyPhen) &
                               !grepl("benign", mutect2_dcs$CLIN_SIG) &
                               !(mutect2_dcs$Variant_Classification %in% c("3'UTR", "5'UTR", "3'Flank", "Intron", "5'Flank", "RNA")), ]
mutect2_dcs <- anti_join(mutect2_dcs, mutect2_all, by = c("Start_Position", "Tumor_Sample_Barcode", "HGVSc"))
mutect2_dcs <- anti_join(mutect2_dcs, mutect2_sscs, by = c("Start_Position", "Tumor_Sample_Barcode", "HGVSc"))

mutect2 <- rbind(mutect2_all, mutect2_sscs, mutect2_dcs)
mutect2 <- mutect2[order(mutect2$Tumor_Sample_Barcode), ]
colnames(mutect2) <- c("hugo", "chr", "start", "end", "variant", "type", "ref", "alt", "dbSNP", "sample", 
                        "HGVSc", "HGVSp", "depth", "ref_count", "alt_count", "SIFT", "PolyPhen", "Clin_Sig")
mutect2$vaf <- mutect2$alt_count/mutect2$depth*100
mutect2 <- mutect2 %>% relocate(vaf, .after = alt_count)
rm(mutect2_all, mutect2_sscs, mutect2_dcs)
write.table(mutect2, file.path(outdir, paste0(project, "_mutect2_", date, ".txt")), sep = "\t", row.names = FALSE)

# Parse mutations (strelka)
strelka_all <- strelka_all[strelka_all$Hugo_Symbol %in% c("BAP1", "GNA11", "GNAQ", "SF3B1", "EIF1AX"), ]
strelka_all <- strelka_all[, colnames(strelka_all) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                                          "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                                          "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG")]
strelka_all <- strelka_all[!grepl("benign", strelka_all$SIFT) &
                             !grepl("benign", strelka_all$PolyPhen) &
                             !grepl("benign", strelka_all$CLIN_SIG) &
                             !(strelka_all$Variant_Classification %in% c("3'UTR", "5'UTR", "3'Flank", "Intron", "5'Flank", "RNA", "Silent", "Splice_Region")), ]

strelka_sscs <- strelka_sscs[strelka_sscs$Hugo_Symbol %in% c("BAP1", "GNA11", "GNAQ", "SF3B1", "EIF1AX"), ]
strelka_sscs <- strelka_sscs[, colnames(strelka_sscs) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                                             "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                                             "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG")]
strelka_sscs <- strelka_sscs[!grepl("benign", strelka_sscs$SIFT) &
                               !grepl("benign", strelka_sscs$PolyPhen) &
                               !grepl("benign", strelka_sscs$CLIN_SIG) &
                               !(strelka_sscs$Variant_Classification %in% c("3'UTR", "5'UTR", "3'Flank", "Intron", "5'Flank", "RNA", "Silent", "Splice_Region")), ]
strelka_sscs <- anti_join(strelka_sscs, strelka_all, by = c("Start_Position", "Tumor_Sample_Barcode", "HGVSc"))

strelka <- rbind(strelka_all, strelka_sscs)
strelka <- strelka[order(strelka$Tumor_Sample_Barcode), ]
colnames(strelka) <- c("hugo", "chr", "start", "end", "variant", "type", "ref", "alt", "dbSNP", "sample", 
                       "HGVSc", "HGVSp", "depth", "ref_count", "alt_count", "SIFT", "PolyPhen", "Clin_Sig")
strelka$vaf <- strelka$alt_count/strelka$depth*100
strelka <- strelka %>% relocate(vaf, .after = alt_count)
rm(strelka_all, strelka_sscs)
write.table(strelka, file.path(outdir, paste0(project, "_strelka_", date, ".txt")), sep = "\t", row.names = FALSE)

# Parse mutations (Pindel)
pindel_all <- pindel_all[pindel_all$Hugo_Symbol %in% c("BAP1"), ]
pindel_all <- pindel_all[, colnames(pindel_all) %in% c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
                                                          "Reference_Allele", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp_Short",
                                                          "t_depth", "t_ref_count", "t_alt_count", "SIFT", "PolyPhen", "CLIN_SIG")]
pindel_all <- pindel_all[!grepl("benign", pindel_all$SIFT) &
                             !grepl("benign", pindel_all$PolyPhen) &
                             !grepl("benign", pindel_all$CLIN_SIG) &
                             !(pindel_all$Variant_Classification %in% c("3'UTR", "5'UTR", "3'Flank", "Intron", "5'Flank", "RNA", "Silent", "Splice_Region")), ]

pindel <- pindel_all[order(pindel_all$Tumor_Sample_Barcode), ]
colnames(pindel) <- c("hugo", "chr", "start", "end", "variant", "type", "ref", "alt", "dbSNP", "sample", 
                       "HGVSc", "HGVSp", "depth", "ref_count", "alt_count", "SIFT", "PolyPhen", "Clin_Sig")
pindel$vaf <- pindel$alt_count/pindel$depth*100
pindel <- pindel %>% relocate(vaf, .after = alt_count)
rm(pindel_all)
