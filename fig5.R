## Combines count and Cnet data in one plot (day 14)
library("dplyr")
library(ggplot2)

# import data
count_bact <- read.csv("data/flow2_day14_bact.csv")
count_pt <- read.csv("data/flow2_day14_pt.csv")
cnet_stat <- read.csv("data/SIP_cnet_day14.csv")

# get statistics from count data
count_bact_stat <- count_bact %>%
  group_by(Treatment, Ring) %>%
  summarize(mean = mean(Abundance), 
            sd = sd(Abundance))

# merbe count and cnet in to df
count_bact_stat$Ring[count_bact_stat$Ring==1] <- 'inner'
count_bact_stat$Ring[count_bact_stat$Ring==2] <- 'outer'
df <- cnet_stat
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
postscript("figures/fig5a.eps", width = 3, height = 3)

ggplot(df, aes(x=cnet_q50, y=count_mean, shape=ring, colour=treatment)) + 
  geom_point() + 
  geom_errorbar(aes(ymax = count_mean+count_sd, ymin = count_mean-count_sd), width=0.02) + 
  geom_errorbarh(aes(xmax = cnet_q75, xmin = cnet_q25), height=0.02) + 
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values=c("#C00000", "#0432FF", "#AB7942", "#000000")) +
  scale_shape_manual(values = c(16,1)) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))