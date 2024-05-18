library("dplyr")
source("utils.R")  # Continued from "fig4_pm_count.R"
detach(package:Rmisc)
detach(package:plyr)


df_fold <- get_df_fold(get_df(), get_df_alg())

df_fold_stat <- df_fold %>%
  group_by(Treatment, Ring) %>%
  summarize(q25 = quantile(rate, probs = 0.25), 
            q50 = quantile(rate, probs = 0.5),
            q75 = quantile(rate, probs = 0.75),
            n = n()
            )


##### PLOT GROWTH RATE
## INNER RING
df_fold_vis <- df_fold[df_fold$Ring==1,]
df_fold_stat_vis <- df_fold_stat[df_fold_stat$Ring==1,]

setEPS()
postscript("figures/fig4_rate_inner.eps", width = 2, height = 2.5)
ggplot() +
  geom_sina(data = df_fold_vis, 
            aes(Treatment, rate, color=Treatment),
            scale = 'width',
            size=1.2,
            maxwidth = 0.6) +
  geom_errorbar(data=df_fold_stat_vis, 
                aes(x=Treatment, ymin=q25, ymax=q75),
                width = 0.15, color='black', size=0.4) + 
  geom_errorbar(data=df_fold_stat_vis, 
                aes(x=Treatment, ymin=q50, ymax=q50),
                width = 0.3, color='black', size=0.8) + 
  scale_color_manual(values=c("#C00000","#0432FF","#AB7942","#000000")) +
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

dev.off()

## stats test inner ring
kruskal.test(rate ~ Treatment, data = df_fold_vis)
pairwise.wilcox.test(
  df_fold_vis$rate, 
  df_fold_vis$Treatment, 
  p.adjust.method = "BH"
)


## OUTER RING
df_fold_vis <- df_fold[df_fold$Ring==2,]
df_fold_stat_vis <- df_fold_stat[df_fold_stat$Ring==2,]

ggplot() +
  geom_sina(data = df_fold_vis, 
            aes(Treatment, rate, color=Treatment),
            scale = 'width',
            size=1.5,
            maxwidth = 0.45, shape=21, fill="white", stroke=0.5) +
  geom_errorbar(data=df_fold_stat_vis, 
                aes(x=Treatment, ymin=q25, ymax=q75),
                width = 0.15, size=0.3, color='black') + 
  geom_errorbar(data=df_fold_stat_vis, 
                aes(x=Treatment, ymin=q50, ymax=q50),
                width = 0.3, size=0.6, color='black') + 
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = 8),
        axis.ticks = element_line(colour = 'black', size=0.2),
        axis.line = element_line(size = 0.2, colour = 'black'),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 8)
  )

ggsave("figures/fig4_rate_outer.pdf", width = 1.5, height = 2.5)

## stats test outer ring
kruskal.test(rate ~ Treatment, data = df_fold_vis)
pairwise.wilcox.test(
  df_fold_vis$rate, 
  df_fold_vis$Treatment, 
  p.adjust.method = "BH"
)


## corrected by algae number
setEPS()
postscript("figures/fig_flow2_outer_corrected.eps")
ggplot() +
  geom_sina(data = df_fold2, aes(Treatment, fold2_log), maxwidth = 0.7) +
  geom_errorbar(data=df_fold2_se, 
                aes(x=Treatment, ymin=fold2_log-se, ymax=fold2_log+se),
                width = 0.3, color='#28A349', size=1) + 
  ylim(-2.5, 5) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()


##### CELL NUMBER
## DAY 5,14 -- inner ring
setEPS()
postscript("figures/fig_count_day5+14_inner.eps", width = 6, height = 5)
boxplot(Abundance~Time+Treatment, 
        data=df[(df$Ring==1)&(df$Treatment!=2),], 
        log = "y",
        # at = c(1:2, 4:5, 7:8),
        col = c('pink','tomato'),
        xlab="Strain in Ring 1", ylab="Bacteria (1/ml)",
        # names=c('Alcani', 'Devosi', 'Marino', 'Empty'),
        main="Inner ring")
dev.off()


## DAY 14 -- outer ring
setEPS()
postscript("figures/fig_count_day14_outer.eps", width = 5, height = 5)
boxplot(Abundance~Treatment, data=df14[df14$Ring==2,],
        xlab="Strain in Ring 1", ylab="Bacteria (1/ml)",
        # names=c('Alcani', 'Devosi', 'Marino', 'Empty'),
        main="Outer ring (Marinobacter)")
dev.off()


## DAY 5 -- outer ring
setEPS()
postscript("figures/fig_count_day5.eps", width = 5, height = 5)
boxplot(Abundance~Treatment, data=df5[df5$Ring==2,],
        xlab="Strain in Ring 1", ylab="Bacteria (1/ml)",
        names=c('Alcani', 'Devosi', 'Marino', 'Empty'),
        main="Outer ring (Marinobacter)")
dev.off()