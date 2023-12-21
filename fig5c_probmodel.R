### Estimates carbon distribution from VB's metabolomics
# 1. import cell count
source("fig4_pm_count.R")
count_df <- get_df()

# 2. get probability matrix
detach(package:Rmisc)
detach(package:plyr)
library("dplyr")

summarize_peak <- function(features_mat){
  list_samples <- rownames(features_mat)
  list_features <- colnames(features_mat)
  samples_df <- data.frame(matrix(nrow=0, ncol=3))
  colnames(samples_df) <- c('sample', 'strain', 'replicate')
  
  for (sample in list_samples){
    str_list <- strsplit(sample, split='-')
    strain <- str_list[[1]][2]
    replicate <- strsplit(str_list[[1]][3], split='_')[[1]][2]
    samples_df[nrow(samples_df)+1, ] <- 
      data.frame(sample = sample, strain = strain, replicate = replicate)
  }
  
  list_strains <- unique(samples_df$strain)
  features_mat_stat <- matrix(nrow=length(list_strains), ncol=length(list_features))
  rownames(features_mat_stat) <- list_strains
  colnames(features_mat_stat) <- list_features
  
  for (strain in list_strains){
    samples_ind <- grep(strain, list_samples)
    features_mat_stat[strain,] <- apply(features_mat[samples_ind,], 2, mean)
  }
  
  return(features_mat_stat)
}


uptake_prob_mat <- data.frame(matrix(nrow=0, ncol=))
colnames(uptake_prob_df) <- c('strain', 'feature', 'P')
list_strains <- unique(features_df_stat$strain)
list_features <- unique(features_df_stat$feature)
peak_control <- features_df_stat[features_df_stat$strain=='abiotic',]
for (strain in list_strains){
  peak_strain <- features_df_stat[features_df_stat$strain==strain,]
}


features_mat <- read.csv("data/lc-ms/peakheights.csv", row.names = 1)
features_mat_mean <- summarize_peak(features_mat)
