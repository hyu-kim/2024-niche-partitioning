source("fig5b.R")

df_summ <- get_total_incorp_df()$df_summ

# plot, average over microplates
setEPS()
postscript("figures/fig5b.eps", width = 1.5, height = 2.5)

df_summ %>% mutate(SD = c_incorp_sd) %>%
  group_by(treatment) %>%
  mutate(SDPos = cumsum(c_incorp_mean)) %>%
  ggplot(aes(x = reorder(treatment, -c_incorp_mean), y = c_incorp_mean, fill = treatment
             # alpha=factor(ring, levels=c("outer","inner"))
             )) +
  geom_bar(stat = "identity", color='black', size=0.2, width=0.7) +
  # geom_text(aes(label = paste(round(c_incorp_mean,3), '±', round(c_incorp_sd,3))), 
  #           position = position_stack(vjust =  0.5)) +
  geom_errorbar(aes(ymin = SDPos-SD, ymax = SDPos+SD), width=0.1, linewidth=0.2, position = "identity", color='black') +
  scale_alpha_discrete(range=c(0.3, 1)) +
  # scale_fill_manual(values=c("#FFDFE1", "#D2DCFB", "#F8DFC4", "#D1D3D4")) +
  scale_fill_manual(values=c("#E06666", "#4061f4", "#AB7942", "#5b5b5b")) +  
  # Alcani, Devosi, Marino, none
  ylab('carbon mass incorporation rate (ng C d-1)') +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        axis.line = element_line(size = 0.2)
        # axis.ticks = element_blank()
  )

dev.off()



## By normalizing with Pt cell number (depreciated)
# df$c_incorp <- df$cnet_q50 * df$count_total * df$mass_c
# df$c_incorp_pt <- df$c_incorp / (df$count_pt * mass[4, 'mass'])
# 
# pdf("figures/fig5b.pdf", width = 2.5, height = 2)
# 
# ggplot(df, aes(x=microplate, 
#                y=c_incorp_pt, 
#                fill=treatment, 
#                alpha=factor(ring, levels=c("outer","inner"))
#                )
#        ) + 
#   geom_bar(stat = 'identity', position = 'stack', color='black', size=0.2) +
#   facet_grid(~factor(treatment, 
#                      levels=c('Devosia', 'Marinobacter', 'Alcanivorax', 'none')
#                      )) +
#   scale_fill_manual(values=c("#E06666", "#4061f4", "#AB7942", "#5b5b5b")) +  # Alcani, Devosi, Marino, none
#   scale_alpha_discrete(range=c(0.3, 1)) +
#   theme(strip.background = element_rect(fill=NA),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         panel.border = element_blank(),
#         legend.position = "None",
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.line = element_line(size = 0.2),
#         strip.text.x = element_blank()
#         # axis.ticks = element_blank()
#         )
# 
# dev.off()