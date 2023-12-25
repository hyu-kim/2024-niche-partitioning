### Estimates carbon distribution from VB's metabolomics
source("fig4_pm_count.R")
detach(package:Rmisc)
detach(package:plyr)
library("dplyr")

get_number_ratio <- function(count_df){
  res_df <- data.frame(matrix(nrow=0, ncol=6))
  colnames(res_df) <- c('Treatment', 'Microplate', 'Direction', 'val1', 'val2', 'ratio')
  list_treatments <- unique(count_df$Treatment)
  
  for (trt in list_treatments){
    # if (trt=='none')
    #   next
    for (mic in seq(6)){
      if (!length(count_df[count_df$Treatment==trt & count_df$Microplate==mic,]))
        next
      for (dir in c(1,3,5)){
        inds <- (count_df$Treatment==trt)&(count_df$Microplate==mic)&(count_df$Direction==dir)
        if (sum(inds) < 2)
          next
        val1 <- count_df$Abundance[inds & count_df$Ring==1] # influencer
        val2 <- count_df$Abundance[inds & count_df$Ring==2] # recipient
        
        res_df[nrow(res_df)+1, ] <- 
          data.frame(Treatment = trt, Microplate = mic, Direction = dir,
                     val1 = val1, val2 = val2, ratio = val2 / val1)
      }
    }
  }
  return(res_df)
}


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


get_uptake_prob <- function(peak_mean){
  uptake_prob <- peak_mean
  uptake_prob[,] <- 0
  
  list_strains <- unique(rownames(peak_mean))
  list_features <- unique(colnames(peak_mean))
  
  peak_control <- peak_mean['abiotic',]
  peak_control <- matrix(rep(peak_control, length(list_strains)), nrow=length(peak_control))
  peak_control <- t(peak_control)
  
  uptake_prob <- (peak_control>peak_mean) * (peak_control-peak_mean) / peak_control
    # reconciliate if this is a correct way to infer probability
    # noting that this is susceptible when there is a false positive in at least one of abiotic control samples
  
  return(uptake_prob)
}


get_distribute_ratio <- function(count_ratio, uptake_prob, exudates_ratio){
  uptake_sums<- rowSums(uptake_prob)
  res_df <- count_ratio[,c('Treatment', 'Microplate', 'Direction', 'ratio')]
  colnames(res_df) <- c('Treatment', 'Microplate', 'Direction', 'count_r')
                        # 'uptake_r', 'exudates_r', 'total_r')
  for (t in unique(res_df$Treatment)){
    if (t == 'none')
      next
    res_df$uptake_r[res_df$Treatment==t] <- uptake_sums['Marinobacter'] / uptake_sums[t]
  }
  
  res_df$exudates_r <- exudates_ratio
  res_df$total_r <- res_df$count_r * res_df$uptake_r * res_df$exudates_r
  
  res_df_stat <- res_df %>%
    group_by(Treatment, Microplate) %>%
    summarize(
      totalmean = mean(total_r),
      totalsd = sd(total_r)
      # geomean = exp(mean(log(total_r))),
      # geosd = exp(sd(log(total_r)))
      )
  
  return(list("df" = res_df, "df_stat" = res_df_stat))
}


get_distributions <- function(peak_loc="data/lc-ms/peakheights.csv", 
                              DOM_loc="data/lc-ms/MicroplateDOM.csv"){
  # 1. import cell count
  count_df <- get_df()
  count_ratio <- get_number_ratio(count_df[count_df$Time==14,])
  
  # 2. get probability matrix
  peak_mat <- read.csv(peak_loc, row.names = 1)
  peak_mean <- summarize_peak(peak_mat)
  uptake_prob <- get_uptake_prob(peak_mean)
  uptake_sums<- rowSums(uptake_prob)
  
  # 3. get DOM concentration
  exudates_df <- read.csv(DOM_loc, row.names = 1)
  exudates_ratio <- exudates_df['14','outer'] / exudates_df['14','inner']
  
  # 4. combine three estimations
  distributions_ratio <- get_distribute_ratio(
    count_ratio, uptake_prob, exudates_ratio)$df_stat
  
  return(distributions_ratio)
}