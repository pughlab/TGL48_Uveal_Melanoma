# extract signature matrix from methylation data matrix with specific signatures
library(tidyverse)
library(arrow)
library(data.table)

#read signatures
# chr-start-end
sig <- read_tsv('tcga_tissue_specific_signature_Eric.tsv')
df <- sig[sig$hyper_cohort == 'TCGA-UVM',c('bin_chr','bin_start','bin_end')]
# Use 'TCGA-LIHC' for liver signatures

df <- df[order(df$bin_chr, df$bin_start),]
df <- as.data.table(df)
setkey(df, bin_chr, bin_start, bin_end)

#read and merge methylation data
setwd('/cluster/projects/pughlab/projects/UVM_cfMeDIP/merged')
data <- read_parquet('uveal.parquet')

#mapping
overlap_result <- foverlaps(
    as.data.table(data),
    as.data.table(df),
    by.x = c('bin_chr', 'bin_start', 'bin_end'),
    nomatch=0
) %>%
  distinct()

# save
print(dim(overlap_result))
write.table(overlap_result, 'uveal.tsv', quote=F, row.names = F, sep='\t')
