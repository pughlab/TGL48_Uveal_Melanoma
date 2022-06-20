# merge medremix resulted data matrix
library(tidyverse)
library(arrow)
library(data.table)

setwd('/cluster/projects/pughlab/projects/UVM_cfMeDIP/cfMeDIP/results/cfmedip_nbglm/')
df <- read.csv('~/Project/cfMeDIP/sample_sheet/samplesheet_uveal_cfmedip.csv')
file_prefix <- c()
for (i in 1:nrow(df)) {
  if (!(df$sample_name[i] %in% file_prefix)) {
    file_prefix <- c(file_prefix, df$sample_name[i])
  }
}

merged_data <- NULL

for (f in file_prefix) {
  p <- paste0(f, '_fit_nbglm.feather')
  data <- read_feather(p)
  data <- data[, c('bin_chr', 'bin_start', 'bin_end', 'methylated_posterior')]
  colnames(data)[ncol(data)] <- f
  if (is.null(merged_data)) {
      merged_data <- data
  } else {
      merged_data <- full_join(merged_data, data, by=c('bin_chr', 'bin_start', 'bin_end'))
  }
}
print(dim(merged_data))

#MeDIP
setwd('/cluster/projects/pughlab/projects/UVM_cfMeDIP/MeDIP/results/cfmedip_nbglm/')

df <- read.csv('~/Project/cfMeDIP/sample_sheet/samplesheet_uveal_medip.csv')
file_prefix <- c()
for (i in 1:nrow(df)) {
  if (!(df$sample_name[i] %in% file_prefix)) {
    file_prefix <- c(file_prefix, df$sample_name[i])
  }
}

for (f in file_prefix) {
  p <- paste0(f, '_fit_nbglm.feather')
  data <- read_feather(p)
  data <- data[, c('bin_chr', 'bin_start', 'bin_end', 'methylated_posterior')]
  colnames(data)[ncol(data)] <- f
  merged_data <- full_join(merged_data, data, by=c('bin_chr', 'bin_start', 'bin_end'))
}

print(dim(merged_data))

HBC <- read_parquet('/cluster/projects/pughlab/projects/CHARM/LFS/Ping_medremix/HBC.parquet')
merged_data <- full_join(merged_data, HBC, by=c('bin_chr', 'bin_start', 'bin_end'))

setwd('/cluster/projects/pughlab/projects/UVM_cfMeDIP/merged')
#save results
write_parquet(merged_data, 'uveal.parquet')
write.table(colnames(merged_data), "sample_names.tsv", quote=F, row.names=F, col.names=F)
