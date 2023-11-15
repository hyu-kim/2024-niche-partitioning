library("phyloseq")
source("r-sparcc/R/sparcc.R")
library(propr) # version 4.2.6

## Import and melt phyloseq object
ps = readRDS("data/hex_16S_phyloseq_enrichment.Rds")
ps_n = transform_sample_counts(ps, 
                               function(x) 100 * x / sum(x)
                               ) # normalize to scale 0-100
df = psmelt(ps_n)  # convert to dataframe

## Subset dataframe to contain relevant ASV only
ind <- (df$OTU=='hex_43') | (df$OTU=='hex_41') | (df$OTU=='hex_60') | 
  (df$OTU=='hex_126') | (df$OTU=='hex_193') | (df$OTU=='hex_5') | 
  (df$OTU=='hex_8')
df_isolates <- df[ind,]

## Extract OTU table and transpose
df <- ps_n@otu_table
df <- t(df)

## run SparCC and downstream analysis
res <- sparcc(df, max.iter = 10)
cor <- res[1]$CORR
asv_vec <- colnames(df)
isolates_vec <- c('hex_43','hex_41','hex_60','hex_126','hex_193','hex_5','hex_8')
