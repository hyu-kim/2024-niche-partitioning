## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
library("ggbreak")

df <- read.csv(file='data/SIP_all.csv')
df_vis <- df[(df$isotope=='c')&(df$distance=='outer'),]
df_vis_stat <- df_vis %>%
  group_by(treatment) %>%
  summarize(q25 = quantile(value, probs = 0.25), 
            q50 = quantile(value, probs = 0.5),
            q75 = quantile(value, probs = 0.75))

# export summarized data
df_stat <- df[df$isotope=='c',] %>%
  group_by(treatment, distance) %>%
  summarize(q25 = quantile(value, probs = 0.25), 
            q50 = quantile(value, probs = 0.5),
            q75 = quantile(value, probs = 0.75))
write.csv(df_stat, 'data/SIP_cnet_summary.csv', row.names=FALSE)


# draw figures
ggplot() +
  geom_sina(data = df_vis, 
            aes(x=treatment, y=value, color=treatment), 
            maxwidth = 0.5, 
            alpha=0.4, 
            size=1) +
  geom_errorbar(data=df_vis_stat, 
                aes(x=treatment, ymin=q25, ymax=q75),
                width = 0.15, color='black', size=0.4) + 
  geom_errorbar(data=df_vis_stat, 
                aes(x=treatment, ymin=q50, ymax=q50),
                width = 0.3, color='black', size=0.8) + 
  labs(x = "Isolate in inner", y = "Incorporation by net",
       title = "Outer ring (Marinobacter)") +
  ylim(NA, 0.24) +
  scale_y_continuous(breaks = append(seq(0, 0.09, 0.03), seq(0.10, 0.25, 0.1))) +
  scale_y_break(c(0.09, 0.10), scales=0.2) +
  scale_color_manual(values=c("#E06666","#4061f4","#7c4202","#5b5b5b")) +  # Alcani, Devosi, Marino, none
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

ggsave("figures/SIP_cnet_day14_outer_break.pdf", width = 4, height = 4)


# statistical test
t.test(df_vis[df_vis$treatment=='Alcanivorax',4], df_vis[df_vis$treatment=='Devosia',4])
t.test(df_vis[df_vis$treatment=='Alcanivorax',4], df_vis[df_vis$treatment=='Marinobacter',4])
t.test(df_vis[df_vis$treatment=='Devosia',4], df_vis[df_vis$treatment=='Marinobacter',4])

t.test(df_vis[df_vis$treatment=='none',4], df_vis[df_vis$treatment=='Alcanivorax',4])
t.test(df_vis[df_vis$treatment=='none',4], df_vis[df_vis$treatment=='Devosia',4])
t.test(df_vis[df_vis$treatment=='none',4], df_vis[df_vis$treatment=='Marinobacter',4])