## Converts count and Cnet into SI indices, compares SI to VB's ECM
# detach(package:Rmisc)
# detach(package:plyr)
library("dplyr")
# library(ggplot2)
# library(scales)
source("fig4_pm_count.R")

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


# extract from ecm matrix
merged_df <- data.frame(
  influencer = c('Alcanivorax', 'Devosia', 'Marinobacter'),
  recipient = c('Marinobacter', 'Marinobacter', 'Marinobacter'),
  ecm = NA,
  mu = NA,
  cnet = NA
)
list_strains_ecm = rownames(ecm_mat)

for (i in 1:nrow(merged_df)){
  r_ecm = grep(merged_df[i, 'recipient'], list_strains_ecm)  # row, recipient
  c_ecm = grep(merged_df[i, 'influencer'], list_strains_ecm)  # col, influencer
  merged_df[i, 'ecm'] <- ecm_mat[r_ecm, c_ecm]
}

# extract count and cnet into merged_df
##blabla
df <- cnet_stat_df
colnames(df) <- c("treatment", "ring", "cnet_q25", "cnet_q50", "cnet_q75")
df$count_mean <- NA
df$count_sd <- NA

for (row in 1:nrow(df)){
  t <- df$treatment[row]
  r <- df$ring[row]
  df[row, 'count_mean'] <- count_bact_stat$mean[count_bact_stat$Treatment==t & count_bact_stat$Ring==r]
  df[row, 'count_sd'] <- count_bact_stat$sd[count_bact_stat$Treatment==t & count_bact_stat$Ring==r]
}


# plot
setEPS()
postscript("figures/fig5a.eps", width = 2, height = 2)

ggplot(df, aes(x=cnet_q50, y=count_mean, shape=ring, colour=treatment)) + 
  geom_point(size=2.5, stroke=0.5) + 
  geom_errorbar(aes(ymax = count_mean+count_sd, ymin = count_mean-count_sd), width=0.02, linewidth=0.2) + 
  geom_errorbarh(aes(xmax = cnet_q75, xmin = cnet_q25), height=0.02, linewidth=0.2) + 
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans='log10', breaks=c(0.1)) +
  scale_y_continuous(trans='log10', breaks=c(10^6, 10^7)) +
  annotation_logticks(short = unit(0.1, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.2, "cm"),
                      size = 0.3) +
  scale_color_manual(values=c("#C00000", "#0432FF", "#AB7942", "#000000")) +  # Alcani, Devosi, Marino, none
  scale_shape_manual(values = c(16,1)) +  # Inner, outer
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_line(colour = "grey80", linewidth=0.1),
        # panel.grid.minor = element_line(colour = "grey80", linewidth=0.2),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_blank()
        )

dev.off()


# plot for each ring
df_sub <- df[df$ring=='outer',]

setEPS()
postscript("figures/fig5a_outer_inset.eps", width = 1.5, height = 1)

ggplot(df_sub, aes(x=cnet_q50, y=count_mean, shape=ring, colour=treatment)) + 
  geom_point(size=2.5, stroke=0.5) + 
  geom_errorbar(aes(ymax = count_mean+count_sd, ymin = count_mean-count_sd), width=0.001, linewidth=0.2) + 
  geom_errorbarh(aes(xmax = cnet_q75, xmin = cnet_q25), height=3e5, linewidth=0.2) + 
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(
    limits = c(0.015, 0.04),
    # trans='log10',
    breaks=seq(0.02, 0.04, 0.005)
    ) +
  # scale_y_continuous(limits = c(2e6, 1e7)
  #                    # trans='log10',
  #                    # breaks=c(10^6, 10^7)
  #                    ) +
  # annotation_logticks(short = unit(0.1, "cm"),
  #                     mid = unit(0.1, "cm"),
  #                     long = unit(0.2, "cm"),
  #                     size = 0.3) +
  scale_color_manual(values=c("#C00000", "#0432FF", "#AB7942", "#000000")) +  # Alcani, Devosi, Marino, none
  scale_shape_manual(values = c(1)) +  # Inner (16), outer (1)
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_line(colour = "grey80", linewidth=0.2),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.2)
        # axis.ticks = element_blank()
  )

dev.off()