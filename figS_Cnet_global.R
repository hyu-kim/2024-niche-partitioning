## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
library("ggbreak")
library(FSA)
source('utils.R')


cnet <- read.csv("data/SIP_cnet_v2.csv")
cnet <- cnet[cnet$Cnet>0,]
cnet$Cnet <- cnet$Cnet*100 # to percent
cnet_info <- read.csv("data/SIP_sample_info.csv")
cnet_append <- append_xnet(cnet, cnet_info)
cnet_append_stat <- cnet_append %>%
  group_by(treatment, ring, microplate) %>%
  summarize(q25 = quantile(Cnet, probs = 0.25),
            q50 = quantile(Cnet, probs = 0.5),
            q75 = quantile(Cnet, probs = 0.75),
            n = n(),
            max = max(Cnet)
  )


# All conditions
ggplot() +
  geom_sina(data = cnet_append,
            aes(x=ring, y=Cnet, color=treatment),
            maxwidth = 0.5,
            alpha=0.3,
            size=0.6) +
  facet_grid(microplate ~ treatment) +
  geom_errorbar(data = cnet_append_stat,
                aes(x=ring, ymin=q50, ymax=q50),
                width = 0.3, size=0.8, color='black') +
  geom_errorbar(data = cnet_append_stat,
                aes(x=ring, ymin=q25, ymax=q75),
                width = 0.15, size=0.4, color='black') +
  geom_text(data = cnet_append_stat,
            aes(x=ring, y=-0.2, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  # scale_y_continuous(breaks = append(seq(0, 9, 3), seq(10, 25, 10))) +
  # scale_y_break(c(9, 10), scales=0.1) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  # Alcani, Devosi, Marino, none
  theme(strip.background = element_rect(fill=NA),
         panel.background = element_rect(fill = "transparent", color = NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.background = element_rect(fill = "transparent", color = NA),
         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
         legend.position = "none",
         text = element_text(colour = 'black', size = 8),
         axis.text = element_text(colour = "black", size = 8),
         axis.ticks = element_line(colour = 'black', size=0.2),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank(),
         axis.line.y.right = element_blank()
  )

ggsave("figures/SIP_Cnet_global.pdf", width = 6, height = 5)


# Outer strain
permil_append_vis <- subset(permil_append, ring=='outer')
permil_app_stat_vis <- subset(permil_app_stat, ring=='outer')

ggplot() +
  geom_sina(data = permil_append_vis,
            aes(x=treatment, y=N_permil, color=treatment),
            maxwidth = 0.5,
            alpha=0.3,
            size=0.6) +
  geom_errorbar(data = permil_app_stat_vis,
                aes(x=treatment, ymin=N_q50, ymax=N_q50),
                width = 0.3, size=0.8, color='black') +
  geom_errorbar(data = permil_app_stat_vis,
                aes(x=treatment, ymin=N_q25, ymax=N_q75),
                width = 0.15, size=0.4, color='black') +
  geom_text(data = permil_app_stat_vis,
            aes(x=treatment, y=-0.2, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  scale_y_continuous(breaks = append(seq(0, 750, 250), seq(1000, 5000, 2000))) +
  scale_y_break(c(750, 751), scales=0.1) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  # Alcani, Devosi, Marino, none
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = 8),
        axis.line = element_line(size = 0.2, colour = 'black'),
        axis.ticks = element_line(colour = 'black', size=0.2),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 8)
  )

ggsave("figures/SIP_N_permil_outer.pdf", width = 3, height = 3)
# ggsave("figures/SIP_N_permil_outer_annot.pdf", width = 3, height = 3)



# # statistical test
# df_outer <- df_vis[df_vis$distance=='outer',]
# df_outer['treatment_p'] <- df_outer['treatment']!='none'
# kruskal.test(value~treatment_p, data=df_outer)
# 
# kruskal.test(value~treatment, data=df_outer)
# dunnTest(value~treatment, data=df_outer, method='holm')
# 
# t.test(df_vis[df_vis$treatment=='Alcanivorax',4], df_vis[df_vis$treatment=='Devosia',4])
# t.test(df_vis[df_vis$treatment=='Alcanivorax',4], df_vis[df_vis$treatment=='Marinobacter',4])
# t.test(df_vis[df_vis$treatment=='Devosia',4], df_vis[df_vis$treatment=='Marinobacter',4])
# 
# t.test(df_vis[df_vis$treatment=='none',4], df_vis[df_vis$treatment=='Alcanivorax',4])
# t.test(df_vis[df_vis$treatment=='none',4], df_vis[df_vis$treatment=='Devosia',4])
# t.test(df_vis[df_vis$treatment=='none',4], df_vis[df_vis$treatment=='Marinobacter',4])