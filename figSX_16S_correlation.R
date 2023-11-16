library("phyloseq")
source("r-sparcc/R/sparcc.R")
library("ggplot2")


## Import and melt phyloseq object
ps = readRDS("data/hex_16S_phyloseq_enrichment.Rds")
ps <- subset_samples(ps, CONDITION=='2')
ps_n = transform_sample_counts(ps, 
                               function(x) 100 * x / sum(x)
                               ) # normalize to scale 0-100
df = psmelt(ps_n)  # convert to dataframe

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

### Plot a correlation heat map
ggplot(data = cor_isol_df, aes(isolate2, isolate1, fill = cor)) +
  geom_tile(color = "white") + 
  geom_text(aes(label=round(cor, digits = 3)), size=3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.5,0.5), space = "Lab", 
                       name="SparCC\nCorrelation") +
  scale_x_discrete(limits=isolates_vec) +
  scale_y_discrete(limits=isolates_vec) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

ggsave("figures/figSX_sparCC_Pt-.pdf", width=6, height=5)
