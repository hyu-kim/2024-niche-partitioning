### Continued from "fig4_pm_count.R"
## PLOT FOLD INCREASE
# Plot by Treatment and Ring
boxplot(fold~Strain, data=df_fold[df_fold$Ring==1,], main="Ring 1", 
        xlab="Strain", ylab="Fold increase")

# before correction, bacteria only
setEPS()
postscript("figures/fig_flow2_outer.eps")
boxplot(fold_log~Treatment, data=df_fold[df_fold$Ring==2,], main="Ring 2 (Marinobacter)", 
        xlab="Strain in Ring 1", ylab="Fold increase",
        names=c('B7WZ', 'EA2', 'None', '3-2')
        # ylim=c(0,20)
)
dev.off()

# after correction with algae number
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

setEPS()
postscript("figures/fig_flow2_outer_corrected_box.eps")
boxplot(fold2_log~Treatment, data=df_fold2,
        xlab="Strain in Ring 1", ylab="Fold increase",
        names=c('B7WZ', 'EA2', 'None', '3-2')
        # ylim=c(0,20)
)
dev.off()


## DAY 5,14 -- inner ring
setEPS()
postscript("figures/fig_count_day5+14_inner.eps", width = 6, height = 5)
boxplot(Abundance~Time+Treatment, 
        data=df[(df$Ring==1)&(df$Treatment!=2),], 
        log = "y",
        at = c(1:2, 4:5, 7:8),
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