library("ggplot2")
library(ggforce)

df <- read.csv(file='data/SIP_roi.csv')
sample_info <- read.csv(file='data/SIP_sample_info.csv')

## add sample details to df using sample_info
for(i in 1:nrow(sample_info)){
  sample_name <- sample_info$sample_name[i]
  df$microplate[df$SAMPLE.NAME==sample_name] <- sample_info$microplate[i]
  df$ring[df$SAMPLE.NAME==sample_name] <- sample_info$ring[i]
  df$strain[df$SAMPLE.NAME==sample_name] <- sample_info$strain[i]
  df$treatment[df$SAMPLE.NAME==sample_name] <- sample_info$treatment[i]
}
df <- df[(df$ring!='inner')|(df$treatment!='none'),]

## summarize df
library("dplyr")
# detach(package:plyr)
df_stat <- df %>%
  group_by(ring, treatment, strain) %>%
  summarize(q25 = quantile(ROIAREA, probs = 0.25), 
            q50 = quantile(ROIAREA, probs = 0.5),
            q75 = quantile(ROIAREA, probs = 0.75),
            mean = mean(ROIAREA)
            )


## plot figures
ggplot() +
  geom_sina(data = df[(df$ring=='inner'),],
            aes(x=treatment, y=ROIAREA, color=treatment),
            maxwidth = 0.5,
            alpha=0.3,
            size=0.8) +
  scale_y_continuous(limits = c(0,5)) +
  geom_errorbar(data = df_stat[df_stat$ring=='inner',],
                aes(x=treatment, ymin=q50, ymax=q50),
                width = 0.3, size=0.8, color='black') +
  geom_errorbar(data=df_stat[df_stat$ring=='inner',],
                aes(x=treatment, ymin=q25, ymax=q75),
                width = 0.15, size=0.4, color='black') +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line = element_line(),
        axis.title = element_blank()
  )

ggsave("figures/fig4d_size_inner.pdf", width = 2.5, height = 2.25)


## statistical test
kruskal.test(ROIAREA ~ ring, data = df[(df$treatment=='Marinobacter'),])
pairwise.wilcox.test(
  df$ROIAREA[(df$ring=='outer')], 
  df$treatment[(df$ring=='outer')], 
  p.adjust.method = "BH"
  )