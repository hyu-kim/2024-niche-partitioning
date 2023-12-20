### Estimates carbon distribution from VB's metabolomics
# import and clean data
source("fig4_pm_count.R")

count_df <- get_df()

features_mat <- read.csv("data/lc-ms/peakheights.csv")