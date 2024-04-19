## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
# detach(package:plyr)
library("ggbreak")
source('fig_s_Xnet_global.R')


nnet <- read.csv("data/SIP_nnet_v2.csv")
nnet <- nnet[nnet$value>0,]
nnet_info <- read.csv("data/SIP_sample_info.csv")
nnet_append <- append_xnet(nnet, nnet_info)
nnet_outer <- nnet_append[(nnet_append$ring=='outer'),]
nnet_outer_stat <- nnet_outer %>%
   group_by(treatment, ring) %>%
   summarize(q25 = quantile(value, probs = 0.25),
             q50 = quantile(value, probs = 0.5),
             q75 = quantile(value, probs = 0.75),
             n = n(),
             max = max(value)
   )


# draw figures
ggplot() +
  geom_sina(data = nnet_outer, 
            aes(x=treatment, y=value, color=treatment), 
            maxwidth = 0.5, 
            alpha=0.4, 
            size=1) +
  geom_errorbar(data=nnet_outer_stat, 
                aes(x=treatment, ymin=q25, ymax=q75),
                width = 0.15, color='black', size=0.4) + 
  geom_errorbar(data=nnet_outer_stat, 
                aes(x=treatment, ymin=q50, ymax=q50),
                width = 0.3, color='black', size=0.8) + 
  geom_text(data=nnet_outer_stat,
            aes(x=treatment, y=max, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  labs(x = "Isolate in inner", y = "value") +
  # ylim(NA, 0.02) +
  scale_y_continuous(breaks = append(seq(0, 0.003, 0.001), seq(0.004, 0.02, 0.008))) +
  scale_y_break(c(0.003, 0.004), scales=0.15) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +  
    # Alcani, Devosi, Marino, none
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

ggsave("figures/SIP_nnet_day14_outer_break.pdf", width = 5, height = 5)


# statistical test
kruskal.test(value ~ treatment, data = nnet_outer)
pairwise.wilcox.test(nnet_outer$value, nnet_outer$treatment, p.adjust.method = "BH")

print(
  cat(
    'Alcanivorax:', (df_vis_stat$q50[1]-df_vis_stat$q50[4]) / df_vis_stat$q50[4], 
    '\n',
    'Devosia:', (df_vis_stat$q50[2]-df_vis_stat$q50[4]) / df_vis_stat$q50[4], 
    '\n',
    'Marinobacter:', (df_vis_stat$q50[3]-df_vis_stat$q50[4]) / df_vis_stat$q50[4],
    '\n'
    )
  )