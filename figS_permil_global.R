## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
library("ggbreak")
library(FSA)
# source('fig5a.R')

append_xnet <- function(xnet, xnet_info){
   xnet_app <- xnet
   xnet_app$microplate <- NaN
   xnet_app$ring <- NaN
   xnet_app$strain <- NaN
   xnet_app$treatment <- NaN
   list_sample_name <- unique(xnet$sample_name)
   for (s in list_sample_name){
      mi <- xnet_info$microplate[xnet_info$sample_name==s]
      ri <- xnet_info$ring[xnet_info$sample_name==s]
      st <- xnet_info$strain[xnet_info$sample_name==s]
      tr <- xnet_info$treatment[xnet_info$sample_name==s]
      
      xnet_app$microplate[xnet_app$sample_name==s] <- mi
      xnet_app$ring[xnet_app$sample_name==s] <- ri
      xnet_app$strain[xnet_app$sample_name==s] <- st
      xnet_app$treatment[xnet_app$sample_name==s] <- tr
   }
   xnet_app <- xnet_app[(xnet_app$treatment!='none' | xnet_app$ring!='inner'),]
   return(xnet_app)
}


cnet <- read.csv("data/SIP_cnet_v2.csv")
cnet <- cnet[cnet$Cnet>0,]
nnet <- read.csv("data/SIP_nnet_v2.csv")
nnet <- nnet[nnet$value>0,]
sample_info <- read.csv("data/SIP_sample_info.csv")

cnet_append <- append_xnet(cnet, sample_info)
nnet_append <- append_xnet(nnet, sample_info)

cnet_append_stat <- cnet_append %>%
  group_by(sample_name, treatment, microplate, ring) %>%
  summarize(q25 = quantile(Cnet, probs = 0.25),
            q50 = quantile(Cnet, probs = 0.5),
            q75 = quantile(Cnet, probs = 0.75),
            n = n(),
            max = max(Cnet)
            # cnet_mean = mean(Cnet),
            # cnet_sd = sd(Cnet),
            # replaced because we look for total incorp rate
  )

nnet_append_stat <- nnet_append %>%
  group_by(treatment, ring) %>%
  summarize(q25 = quantile(value, probs = 0.25),
            q50 = quantile(value, probs = 0.5),
            q75 = quantile(value, probs = 0.75),
            n = n(),
            max = max(value)
  )


# draw figures for Cnet
ggplot() +
  geom_sina(data = cnet_append,
            aes(x=ring, y=value, color=treatment),
            maxwidth = 0.8,
            alpha=0.4,
            size=0.3) +
  facet_grid(cols = vars(treatment), rows = vars(microplate)) +
  geom_errorbar(data=cnet_append_stat,
                aes(x=ring, ymin=cnet_q25, ymax=cnet_q75),
                width = 0.15, color='black', size=0.4) +
  geom_errorbar(data=cnet_append_stat,
                aes(x=ring, ymin=cnet_q50, ymax=cnet_q50),
                width = 0.3, color='black', size=0.8) +
  geom_text(data=cnet_append_stat,
            aes(x=ring, y=cnet_max, label=paste('n = ', cnet_n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  labs(y = expression(C[net]), x = "Location") +
  # ylim(NA, 0.24) +
  # scale_y_continuous(breaks = append(seq(0, 0.09, 0.03), seq(0.10, 0.25, 0.1))) +
  # scale_y_break(c(0.09, 0.10), scales=0.2) +
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

# ggsave("figures/SIP_cnet_global_v3.pdf", width = 8, height = 5)


# figures for Nnet
ggplot() +
  geom_sina(data = nnet_append,
            aes(x=ring, y=100*value, color=treatment),
            maxwidth = 0.8,
            alpha=0.3,
            size=0.5) +
  facet_grid(cols = vars(treatment)) +
  geom_errorbar(data=nnet_append_stat,
                aes(x=ring, ymin=100*q25, ymax=100*q75),
                width = 0.15, color='black', size=0.4) +
  geom_errorbar(data=nnet_append_stat,
                aes(x=ring, ymin=100*q50, ymax=100*q50),
                width = 0.3, color='black', size=0.8) +
  geom_text(data=nnet_append_stat,
            aes(x=ring, y=100*max, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  labs(y = expression(N[net]), x = "Location") +
  scale_y_continuous(breaks = append(seq(0, 0.3, 0.1), seq(0.3, 2, 0.5))) +
  scale_y_break(c(0.3, 0.32), scales=0.15) +
  # ylim(NA, 0.01) +
  # scale_y_continuous(breaks = append(seq(0, 0.09, 0.03), seq(0.10, 0.25, 0.1))) +
  # scale_y_break(c(0.09, 0.10), scales=0.2) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  # Alcani, Devosi, Marino, none
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.text = element_text(colour = "black")
  )

ggsave("figures/SIP_nnet_global.pdf", width = 6, height = 4)


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