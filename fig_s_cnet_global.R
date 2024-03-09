## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
library("ggbreak")
library(FSA)
source('fig5a.R')

append_cnet <- function(cnet, cnet_info){
   cnet_app <- cnet
   cnet_app$microplate <- NaN
   cnet_app$ring <- NaN
   cnet_app$strain <- NaN
   cnet_app$treatment <- NaN
   list_sample_name <- unique(cnet$sample_name)
   for (s in list_sample_name){
      mi <- cnet_info$microplate[cnet_info$sample_name==s]
      ri <- cnet_info$ring[cnet_info$sample_name==s]
      st <- cnet_info$strain[cnet_info$sample_name==s]
      tr <- cnet_info$treatment[cnet_info$sample_name==s]
      
      cnet_app$microplate[cnet_app$sample_name==s] <- mi
      cnet_app$ring[cnet_app$sample_name==s] <- ri
      cnet_app$strain[cnet_app$sample_name==s] <- st
      cnet_app$treatment[cnet_app$sample_name==s] <- tr
   }
   cnet_app <- cnet_app[(cnet_app$treatment!='none' | cnet_app$ring!='inner'),]
   return(cnet_app)
}

cnet_append <- append_cnet(cnet, cnet_info)

cnet_append_stat <- cnet_append %>%
   group_by(sample_name, treatment, microplate, ring) %>%
   summarize(cnet_q25 = quantile(Cnet, probs = 0.25),
             cnet_q50 = quantile(Cnet, probs = 0.5),
             cnet_q75 = quantile(Cnet, probs = 0.75),
             cnet_n = n(),
             cnet_max = max(Cnet)
             # cnet_mean = mean(Cnet),
             # cnet_sd = sd(Cnet),
             # replaced because we look for total incorp rate
   )


# draw figures
ggplot() +
  geom_sina(data = cnet_append, 
            aes(x=ring, y=Cnet, color=treatment), 
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

ggsave("figures/SIP_cnet_global_v3.pdf", width = 8, height = 5)


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