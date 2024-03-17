## Combines count and Cnet data in one plot (day 14)
library("dplyr")
library(ggplot2)
library(scales)

summarize_cnet <- function(cnet, cnet_info){
  # summarize cnet
  cnet$microplate <- NA
  cnet$ring <- NA
  cnet$treatment <- NA
  cnet$strain <- NA
  for (s in unique(cnet_info$sample_name)){
    m <- cnet_info$microplate[cnet_info$sample_name==s]
    r <- cnet_info$ring[cnet_info$sample_name==s]
    t <- cnet_info$treatment[cnet_info$sample_name==s]
    st <- cnet_info$strain[cnet_info$sample_name==s]
    
    cnet$microplate[cnet$sample_name==s] <- m
    cnet$ring[cnet$sample_name==s] <- r
    cnet$treatment[cnet$sample_name==s] <- t
    cnet$strain[cnet$sample_name==s] <- st
  }
  
  cnet_summ <- cnet %>%
    group_by(treatment, microplate, ring) %>%
    summarize(cnet_q25 = quantile(Cnet, probs = 0.25),
              cnet_q50 = quantile(Cnet, probs = 0.5),
              cnet_q75 = quantile(Cnet, probs = 0.75),
              cnet_n = n()
              # cnet_mean = mean(Cnet),
              # cnet_sd = sd(Cnet),
              # replaced because we look for total incorp rate
    )
  return(list('df' = cnet, 'df_summ' = cnet_summ))
}


merge_count_cnet <- function(count_bact_stat, cnet_stat){
  count_bact_stat$Ring[count_bact_stat$Ring==1] <- 'inner'
  count_bact_stat$Ring[count_bact_stat$Ring==2] <- 'outer'
  df <- cnet_stat
  colnames(df) <- c("treatment", "microplate", "ring", "cnet_q25", "cnet_q50", 
                    "cnet_q75", "n_cnet")
  df$count_mean <- NA
  df$count_sd <- NA
  
  for (row in 1:nrow(df)){
    t <- df$treatment[row]
    r <- df$ring[row]
    m <- df$microplate[row]
    df[row, 'count_mean'] <- 
      count_bact_stat$mean[count_bact_stat$Treatment==t 
                           & count_bact_stat$Ring==r
                           & count_bact_stat$Microplate==m]
    df[row, 'count_sd'] <- 
      count_bact_stat$sd[count_bact_stat$Treatment==t 
                         & count_bact_stat$Ring==r
                         & count_bact_stat$Microplate==m]
  }
  
  df <- df[(df$treatment!='none' | df$ring!='inner'),]
  
  return(df)
}

# import and clean data
count_bact <- read.csv("data/flow2_day14_bact.csv")
count_pt <- read.csv("data/flow2_day14_pt.csv")
cnet <- read.csv("data/SIP_cnet_v2.csv")
cnet <- cnet[cnet$Cnet>0,]
cnet_info <- read.csv("data/SIP_sample_info.csv")

# summarize cnet
cnet_stat <- summarize_cnet(cnet, cnet_info)$df_summ

# get statistics from count data
count_bact_stat <- count_bact %>%
  group_by(Treatment, Microplate, Ring) %>%
  summarize(mean = mean(Abundance), 
            sd = sd(Abundance))

# merge count and cnet in to df
merged_stat <- merge_count_cnet(count_bact_stat, cnet_stat)


# plot
setEPS()
# postscript("figures/fig5a_lagend.eps", width = 6, height = 6)
postscript("figures/fig5a_v2.eps", width = 1.6, height = 1.6)

ggplot(merged_stat[merged_stat$ring=='inner',], aes(x=cnet_q50, y=count_mean, colour=treatment)) + 
  geom_point(aes(size=n_cnet, fill=treatment, shape=ring), stroke=0.5) + 
  geom_errorbar(aes(ymax = count_mean+count_sd, ymin = count_mean-count_sd), width=0, linewidth=0.1) + 
  geom_errorbarh(aes(xmax = cnet_q75, xmin = cnet_q25), height=0, linewidth=0.1) + 
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans='log10', breaks=c(0.1)) +
  scale_y_continuous(trans='log10', breaks=c(10^6, 10^7)) +
  annotation_logticks(short = unit(0.1, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.2, "cm"),
                      size = 0.1) +
  scale_size(range = c(1, 2.5)) +
  scale_color_manual(values=c("#C00000", "#0432FF", "#AB7942", "#000000")) +
  scale_fill_manual(values=c("#FFDFE1", "#D2DCFB", "#F8DFC4", "#ffffff")) +
    # Alcani, Devosi, Marino, none
  scale_shape_manual(values = c(16,21)) +  # Inner, outer
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_line(colour = "grey80", linewidth=0.2),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.15),
        axis.ticks = element_blank()
        )

dev.off()


# plot outer ring
merged_stat_outer <- merged_stat[merged_stat$ring=='outer',]

setEPS()
postscript("figures/fig5b_v3.eps", width = 1.6, height = 1.6)

ggplot(merged_stat_outer, aes(x=cnet_q50, y=count_mean, shape=ring, colour=treatment)) + 
  geom_point(aes(size=n_cnet, fill=treatment, shape=ring), stroke=0.5) + 
  geom_errorbar(aes(ymax = count_mean+count_sd, ymin = count_mean-count_sd), width=0, linewidth=0.1) +
  geom_errorbarh(aes(xmax = cnet_q75, xmin = cnet_q25), height=0, linewidth=0.1) +
  scale_y_continuous(limits = c(0, 12.5e6), breaks=seq(0, 12e6, 4e6)) +
  # scale_x_continuous(limits = c(0.015, 0.05), breaks=seq(0.02, 0.05, 0.01)) +
  scale_size(range = c(1, 2.5)) +
  scale_color_manual(values=c("#C00000", "#0432FF", "#AB7942", "#000000")) +
  scale_fill_manual(values=c("#FFDFE1", "#D2DCFB", "#F8DFC4", "#ffffff")) +
  # scale_fill_manual(values=c("#FFDFE1", "#D2DCFB", "#F8DFC4", "#D1D3D4")) +
    # Alcani, Devosi, Marino, none
  scale_shape_manual(values = c(21)) +  # Inner, outer
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_line(colour = "grey80", linewidth=0.2),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.15),
        axis.ticks = element_line(size = 0.15)
  )

dev.off()