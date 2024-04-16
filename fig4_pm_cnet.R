## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
# detach(package:plyr)
library("ggbreak")
source('fig_s_Xnet_global.R')


cnet <- read.csv("data/SIP_cnet_v2.csv")
cnet <- cnet[cnet$Cnet>0,]
cnet_info <- read.csv("data/SIP_sample_info.csv")
cnet_append <- append_cnet(cnet, cnet_info)
cnet_outer <- cnet_append[(cnet_append$ring=='outer'),]
cnet_outer_stat <- cnet_outer %>%
   group_by(treatment, ring) %>%
   summarize(q25 = quantile(Cnet, probs = 0.25),
             q50 = quantile(Cnet, probs = 0.5),
             q75 = quantile(Cnet, probs = 0.75),
             n = n(),
             max = max(Cnet)
   )


# draw figures
ggplot() +
  geom_sina(data = cnet_outer, 
            aes(x=treatment, y=Cnet, color=treatment), 
            maxwidth = 0.5, 
            alpha=0.3, 
            size=1) +
  geom_errorbar(data=cnet_outer_stat, 
                aes(x=treatment, ymin=q25, ymax=q75),
                width = 0.15, color='black', size=0.4) + 
  geom_errorbar(data=cnet_outer_stat, 
                aes(x=treatment, ymin=q50, ymax=q50),
                width = 0.3, color='black', size=0.8) + 
  labs(x = "Isolate in inner", y = "Cnet") +
  ylim(NA, 0.24) +
  scale_y_continuous(breaks = append(seq(0, 0.09, 0.03), seq(0.10, 0.25, 0.1))) +
  scale_y_break(c(0.09, 0.10), scales=0.2) +
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

ggsave("figures/SIP_cnet_day14_outer_break_v2.pdf", width = 3, height = 4)
ggsave("figures/SIP_cnet_day14_inner_break.pdf", width = 3, height = 4)


# statistical test
kruskal.test(Cnet ~ treatment, data = cnet_outer)
pairwise.wilcox.test(cnet_outer$Cnet, cnet_outer$treatment, p.adjust.method = "BH")

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