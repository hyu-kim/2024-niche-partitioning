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


permil_df <- read.csv("data/SIP_permil_v2.csv", check.names = FALSE)
sample_info <- read.csv("data/SIP_sample_info.csv")
sample_info$microplate <- as.character(sample_info$microplate)

permil_append <- append_xnet(permil_df, sample_info)
cnet_append <- append_xnet(cnet, sample_info)
nnet_append <- append_xnet(nnet, sample_info)

permil_app_stat <- permil_append %>%
  group_by(sample_name, treatment, microplate, ring) %>%
  summarize(C_q25 = quantile(C_permil, probs = 0.25),
            C_q50 = quantile(C_permil, probs = 0.5),
            C_q75 = quantile(C_permil, probs = 0.75),
            N_q25 = quantile(N_permil, probs = 0.25),
            N_q50 = quantile(N_permil, probs = 0.5),
            N_q75 = quantile(N_permil, probs = 0.75),
            n = n()
            # cnet_mean = mean(Cnet),
            # cnet_sd = sd(Cnet),
            # replaced because we look for total incorp rate
  )


# figures for inner strains
plot_scatter <- function(df=permil_append, loc='inner', save=TRUE){
  df_vis <- subset(df, ring==loc)
  p <- ggplot() +
    geom_point(data = df_vis,
               aes(x=N_permil, y=C_permil, color=treatment),
               alpha=0.4,
               size=1.5) +
    # scale_y_continuous(trans = 'log10') +
    # scale_x_continuous(trans = 'log10') +
    
  # facet_grid(cols = vars(treatment), rows = vars(microplate)) +
  # geom_errorbar(data=cnet_append_stat,
  #               aes(x=ring, ymin=cnet_q25, ymax=cnet_q75),
  #               width = 0.15, color='black', size=0.4) +
  # geom_errorbar(data=cnet_append_stat,
  #               aes(x=ring, ymin=cnet_q50, ymax=cnet_q50),
  #               width = 0.3, color='black', size=0.8) +
  # geom_text(data=cnet_append_stat,
  #           aes(x=ring, y=cnet_max, label=paste('n = ', cnet_n, sep='')),
  #           position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
    labs(x = '15N permil', y = "13C permil", title=loc) +
    scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
    # Alcani, Devosi, Marino, none
    theme(strip.background = element_rect(fill=NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          legend.position = "bottom",
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.line.y.right = element_blank()
    )
  if(loc=='outer'){
    p + ylim(NA, 1200)+ xlim(NA, 750)
  }
    
  if(save){
    ggsave(paste("figures/SIP_permil_",loc,".pdf",sep=''), width = 6, height = 6)
  }
  return()
}

plot_scatter(permil_append, 'inner', TRUE)
plot_scatter(permil_append, 'outer', TRUE)

permil_append_vis <- subset(permil_append, ring=='inner')
ggplot() +
  geom_point(data = permil_append_vis,
            aes(x=N_permil, y=C_permil, color=treatment),
            alpha=0.5,
            size=1) +
  # facet_grid(cols = vars(treatment), rows = vars(microplate)) +
  # geom_errorbar(data=cnet_append_stat,
  #               aes(x=ring, ymin=cnet_q25, ymax=cnet_q75),
  #               width = 0.15, color='black', size=0.4) +
  # geom_errorbar(data=cnet_append_stat,
  #               aes(x=ring, ymin=cnet_q50, ymax=cnet_q50),
  #               width = 0.3, color='black', size=0.8) +
  # geom_text(data=cnet_append_stat,
  #           aes(x=ring, y=cnet_max, label=paste('n = ', cnet_n, sep='')),
  #           position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  labs(y = '15N permil', x = "13c permil", title='inner') +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
    # Alcani, Devosi, Marino, none
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "bottom",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank()
        )

ggsave("figures/SIP_permil_inner.pdf", width = 4, height = 6)


# figures for N-permil of Marinobacter
permil_app_stat2 <- permil_append %>%
  group_by(treatment, ring, strain) %>%
  summarize(C_q25 = quantile(C_permil, probs = 0.25),
            C_q50 = quantile(C_permil, probs = 0.5),
            C_q75 = quantile(C_permil, probs = 0.75),
            C_min = min(C_permil),
            N_q25 = quantile(N_permil, probs = 0.25),
            N_q50 = quantile(N_permil, probs = 0.5),
            N_q75 = quantile(N_permil, probs = 0.75),
            N_min = min(N_permil),
            n = n()
  )

permil_append_vis <- subset(permil_append, strain=='Marinobacter')
permil_app_stat2_vis <- subset(permil_app_stat2, strain=='Marinobacter')


ggplot() +
  geom_sina(data = permil_append_vis,
            aes(x=ring, y=N_permil, color=treatment),
            maxwidth = 0.8,
            alpha=0.5,
            size=0.4) +
  facet_grid(cols = vars(treatment)) +
  geom_errorbar(data=permil_app_stat2_vis,
                aes(x=ring, ymin=N_q25, ymax=N_q75),
                width = 0.15, color='black', size=0.4) +
  geom_errorbar(data=permil_app_stat2_vis,
                aes(x=ring, ymin=N_q50, ymax=N_q50),
                width = 0.3, color='black', size=0.8) +
  geom_text(data=permil_app_stat2_vis,
            aes(x=ring, y=N_min, label=paste('n = ', n, sep='')),
            position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  scale_y_continuous(breaks = append(seq(0, 750, 250), seq(1000, 5000, 1000))) +
  scale_y_break(c(750, 751), scales=0.2) +
  # scale_y_continuous(breaks = seq(0, 5000, 250)) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  # Alcani, Devosi, Marino, none
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        legend.position = "none",
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = 'black', size=0.2),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        text = element_text(size = 8)
  )

ggsave("figures/SIP_N_permil_Marino.pdf", width = 5, height = 4)


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