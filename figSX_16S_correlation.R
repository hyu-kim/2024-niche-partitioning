library("phyloseq")
source("r-sparcc/R/sparcc.R")
library(propr) # version 4.2.6

## Import and melt phyloseq object
ps = readRDS("data/hex_16S_phyloseq_enrichment.Rds")
ps <- subset_samples(ps, CONDITION=='1')
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

## run SparCC and define gene of interest
res <- sparcc(df, max.iter = 10)
cor <- res[1]$CORR
asv_vec <- colnames(df)
isolates_asv_vec <- c('hex_43','hex_41','hex_60','hex_126','hex_193','hex_5','hex_8')
isolates_vec <- c('Marinobacter 3-2', 'Alcanivorax EA2', 'Thalassospira 13M1', 
                'Devosia EAB7WZ', 'Roseibium 13C1', 'Oceanicaulis 13A',
                'Algoriphagus ARW1R1')

## Create a dataframe from correlation
cor_isol_df <- data.frame(matrix(ncol=5, nrow=0))
colnames(cor_isol_df) <- c('hex1', 'hex2', 'isolate1', 'isolate2', 'cor')

for (i in 1:(length(isolates_vec)-1)){
  for (j in (i+1):length(isolates_vec)){
    hex1 <- isolates_asv_vec[i]
    hex2 <- isolates_asv_vec[j]
    isolate1 <- isolates_vec[i]
    isolate2 <- isolates_vec[j]
    ind1 <- which(asv_vec==hex1)
    ind2 <- which(asv_vec==hex2)
    val <- cor[ind1, ind2]
    
    df_app <- data.frame(hex1=isolates_asv_vec[i], hex2=isolates_asv_vec[j], 
                         isolate1=isolates_vec[i], isolate2=isolates_vec[j],
                         cor=val)
    cor_isol_df <- rbind(cor_isol_df, df_app)
  }
}



