## Converts count and Cnet into SI indices, compares SI to VB's ECM
source("fig4_pm_count.R")
detach(package:Rmisc)
detach(package:plyr)
library("dplyr")
# library(ggplot2)
library(scales)

# import data
count_df <- get_df()
count_df_alg <- get_df_alg()
count_df_fold <- get_df_fold(count_df, count_df_alg)
cnet_df <- read.csv("data/SIP_cnet_summary.csv")
ecm_mat <- read.csv("data/lc-ms/ExpectedCompetition.csv", row.names = 1)

# summarize growth rate
count_df_fold_stat <- count_df_fold %>%
  group_by(Treatment, Ring) %>%
  summarize(q25 = quantile(rate, probs = 0.25), 
            q50 = quantile(rate, probs = 0.5),
            q75 = quantile(rate, probs = 0.75))
count_df_fold_stat$Ring[count_df_fold_stat$Ring==1] <- 'inner'
count_df_fold_stat$Ring[count_df_fold_stat$Ring==2] <- 'outer'


# merge ECM, SI from PM
merged_df <- data.frame(
  influencer = rep(c('Alcanivorax', 'Devosia', 'Marinobacter', 'none'), 3),
  recipient = 'Marinobacter',
  measure = c(rep('ecm', 4), rep('pm_mu', 4), rep('pm_cnet', 4)),
  val = 0
)
list_strains_ecm = rownames(ecm_mat)

for (i in 1:nrow(merged_df)){
  influencer_i <- merged_df[i, 'influencer']
  
  if (influencer_i=='none'){
    next
  }
  
  if (merged_df[i, 'measure'] == 'ecm'){
    r_ecm = grep('Marinobacter', list_strains_ecm)  # row, recipient
    c_ecm = grep(influencer_i, list_strains_ecm)  # col, influencer
    merged_df[i, 'val'] <- ecm_mat[r_ecm, c_ecm]
  } else if (merged_df[i, 'measure'] == 'pm_mu') {
    mu_i <- count_df_fold_stat$q50[
      count_df_fold_stat$Treatment==influencer_i & 
        count_df_fold_stat$Ring=='outer'
      ]
    mu_self <- count_df_fold_stat$q50[
      count_df_fold_stat$Treatment=='Marinobacter' & 
        count_df_fold_stat$Ring=='outer'
      ]
    mu_0 <- count_df_fold_stat$q50[count_df_fold_stat$Treatment=='none']
    merged_df[i, 'val'] <- - (mu_i-mu_0) / (mu_self-mu_0)
  } else {
    cnet_i <- cnet_df$q50[cnet_df$treatment==influencer_i & 
                            cnet_df$distance=='outer']
    cnet_self <- cnet_df$q50[cnet_df$treatment=='Marinobacter' &
                               cnet_df$distance=='outer']
    cnet_0 <- cnet_df$q50[cnet_df$treatment=='none']
    merged_df[i, 'val'] <- - (cnet_i - cnet_0) / (cnet_self - cnet_0)
  }
}


# plot by heat
setEPS()
postscript("figures/fig6_summary.eps", width = 3.5, height = 3)

ggplot(merged_df, aes(x=influencer, y=-val, color=measure)) + 
  geom_bar(width = 0.8, stat='identity', position=position_dodge(width=.9), fill='grey95') + 
  # scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_line(colour = "grey80", linewidth=0.1),
        # panel.grid.minor = element_line(colour = "grey80", linewidth=0.2),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(3, 'mm'),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        text = element_text(size = 7),
        axis.line = element_line(size = 0.25),
        )

dev.off()
