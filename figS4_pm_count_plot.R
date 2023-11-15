detach(package:Rmisc)
detach(package:plyr)
library("dplyr")

### Continued from "fig4_pm_count.R"
# bacteria growth rate
df_fold$rate <- df_fold$fold_log / 9
df_fold_stat <- df_fold %>%
  group_by(Treatment, Ring) %>%
  summarize(q25 = quantile(rate, probs = 0.25), 
            q50 = quantile(rate, probs = 0.5),
            q75 = quantile(rate, probs = 0.75))

df_fold$Ring[df_fold$Ring==1] <- 'inner'
df_fold$Ring[df_fold$Ring==2] <- 'outer'


##### PLOT GROWTH RATE
## ALL RINGS
df_fold_vis <- df_fold[df_fold$Ring==1,]
df_fold_stat_vis <- df_fold_stat[df_fold_stat$Ring==1,]

ggplot() +
  # geom_errorbar(data=df_vis_stat, 
  #               aes(x=distance, ymin=q25, ymax=q75),
  #               width = 0.15, color='black', size=0.4) + 
  # geom_errorbar(data=df_vis_stat, 
  #               aes(x=distance, ymin=q50, ymax=q50),
  #               width = 0.3, color='black', size=0.8) + 
  geom_sina(data = df_fold, 
            aes(x=Ring, y=rate, color=Treatment),
            scale = 'width',
            size=1.2,
            maxwidth = 0.6) +
  facet_grid(cols = vars(Treatment)) +
  geom_errorbar(data=df_fold_stat, 
                aes(x=Ring, ymin=q25, ymax=q75),
                width = 0.15, color='black', size=0.4) + 
  geom_errorbar(data=df_fold_stat, 
                aes(x=Ring, ymin=q50, ymax=q50),
                width = 0.3, color='black', size=0.8) + 
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  labs(y = "Growth rate (d-1)", x = "Location") +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank()
  )

ggsave("figures/figS4_rate_global.pdf", width = 6, height = 4)