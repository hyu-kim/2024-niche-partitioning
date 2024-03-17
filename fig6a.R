## Scales SI, ECM, MRO, and total incorp rate to 0-1
source("fig4_pm_count.R")
detach(package:Rmisc)
detach(package:plyr)
library("dplyr")
# library(ggplot2)
library(scales)
source("fig5c.R")


merge_indices <- function(ecm_mat, mro_mat, incorp_df){
  merged_df <- data.frame(
    influencer = rep(c('Alcanivorax', 'Devosia', 'Marinobacter', 'none'), 3),
    recipient = 'Marinobacter',
    measure = c(rep('ecm', 4), rep('mro', 4), rep('c_incorp', 4)),
    val = 0
  )
  merged_df$measure = factor(merged_df$measure, levels = c("ecm", "mro", "c_incorp"), ordered = TRUE)
  list_strains_ecm = rownames(ecm_mat)
  list_strains_mro = colnames(mro_mat)
  
  for (i in 1:nrow(merged_df)){
    influencer_i <- merged_df[i, 'influencer']
    
    if (influencer_i=='none'){
      next
    }
    
    if (merged_df[i, 'measure'] == 'ecm'){
      r_ecm = grep('Marinobacter', list_strains_ecm)  # row, recipient
      c_ecm = grep(influencer_i, list_strains_ecm)  # col, influencer
      merged_df[i, 'val'] <- ecm_mat[r_ecm, c_ecm]
      print('ecm done')
    } else if (merged_df[i, 'measure'] == 'mro'){
      r_ecm = grep('3.2', list_strains_mro)  # row, recipient
      c_ecm = grep(influencer_i, list_strains_mro)  # col, influencer
      merged_df[i, 'val'] <- - mro_mat[r_ecm, c_ecm]
      print('mro done')
    } else if (merged_df[i, 'measure'] == 'c_incorp') {
      incorp_i <- incorp_df$c_incorp_mean[
        incorp_df$treatment==influencer_i & incorp_df$ring=='outer'
      ]
      incorp_self <- incorp_df$c_incorp_mean[
        incorp_df$treatment=='Marinobacter' & incorp_df$ring=='outer'
      ]
      incorp_0 <- incorp_df$c_incorp_mean[incorp_df$treatment=='none']
      merged_df[i, 'val'] <- - (incorp_i-incorp_0) / (incorp_self-incorp_0)
      print('incorp done')
    }
  }
  return(merged_df)
}


# import data
incorp_df <- get_total_incorp_df()$df_summ
# count_df <- get_df()
# count_df_alg <- get_df_alg()
# count_df_fold <- get_df_fold(count_df, count_df_alg)
ecm_mat <- read.csv("data/lc-ms/ExpectedCompetition.csv", row.names = 1)
mro_mat <- read.csv("data/BiologMe/MRO_table.txt", row.names=1) # influencer, each column
merged_df <- merge_indices(ecm_mat, mro_mat, incorp_df)

# # summarize growth rate
# count_df_fold_stat <- count_df_fold %>%
#   group_by(Treatment, Ring) %>%
#   summarize(q25 = quantile(rate, probs = 0.25), 
#             q50 = quantile(rate, probs = 0.5),
#             q75 = quantile(rate, probs = 0.75))
# count_df_fold_stat$Ring[count_df_fold_stat$Ring==1] <- 'inner'
# count_df_fold_stat$Ring[count_df_fold_stat$Ring==2] <- 'outer'




# plot by heat
merged_df_vis <- merged_df[merged_df$influencer=='Alcanivorax' | merged_df$influencer=='Devosia',]

ggplot(merged_df_vis, aes(x=factor(measure, level=unique(measure)), y=-val, color=influencer)) + 
  geom_bar(width = 0.8, stat='identity', position=position_dodge(width=.9), fill='grey95') + 
  ylab('Competition index [AU]') +
  # scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "right",
        legend.key.size = unit(3, 'mm'),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        text = element_text(size = 7),
        axis.line = element_line(size = 0.25),
        )

ggsave("figures/fig6a_v2_draft.pdf", width = 3.5, height = 3)
