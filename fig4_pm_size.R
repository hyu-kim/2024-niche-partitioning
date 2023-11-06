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

## summarize df
library("dplyr")
# detach(package:plyr)
df_stat <- df %>%
  group_by(ring, treatment) %>%
  summarize(q25 = quantile(ROIAREA, probs = 0.25), 
            q50 = quantile(ROIAREA, probs = 0.5),
            q75 = quantile(ROIAREA, probs = 0.75),
            mean = mean(ROIAREA)
            )


## plot figures
ggplot() +
  geom_sina(data = df[(df$ring=='outer'),], 
            aes(x=treatment, y=ROIAREA, color=treatment), 
            maxwidth = 0.5, 
            alpha=0.4, 
            size=1) +
  geom_errorbar(data=df_stat[df_stat$ring=='outer',], 
                aes(x=treatment, ymin=q25, ymax=q75),
                width = 0.15, color='black', size=0.4) + 
  geom_errorbar(data=df_stat[df_stat$ring=='outer',], 
                aes(x=treatment, ymin=q50, ymax=q50),
                width = 0.3, color='black', size=0.8) + 
  labs(x = "Isolate in inner", y = "Incorporation by net",
       title = "Outer ring (Marinobacter)")