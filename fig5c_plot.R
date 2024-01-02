source('fig5c_probmodel.R')

total_incorp_ratio <- get_ratio_df()

# plot results

setEPS()
postscript("figures/fig5c.eps", width = 1.7, height = 1.7)

ggplot(total_incorp_ratio, aes(x=pred_mean, y=incorp_ratio, colour=Treatment)) + 

  geom_point(size=2, stroke=0.5) + 
  geom_errorbarh(aes(xmax = pred_mean+pred_sd, xmin = pred_mean-pred_sd), height=0.02, linewidth=0.2) + 
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  annotation_logticks(short = unit(0.1, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.2, "cm"),
                      size = 0.15) +
  scale_color_manual(values=c("#C00000", "#0432FF", "#AB7942", "#000000")) +  # Alcani, Devosi, Marino, none
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_line(colour = "grey80", linewidth=0.2),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2)
  )

dev.off()