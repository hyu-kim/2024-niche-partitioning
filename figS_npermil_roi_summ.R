library("dplyr")
library("ggplot2")
source('utils.R')

permil_df <- read.csv("data/SIP_permil_v2.csv", check.names = FALSE)
sample_info <- read.csv("data/SIP_sample_info.csv")
sample_info$microplate <- as.character(sample_info$microplate)

permil_append <- append_xnet(permil_df, sample_info)

permil_app_stat <- permil_append %>%
  group_by(strain, treatment, ring, microplate) %>%  # better to exclude microplate for visibility..
  summarize(N_q25 = quantile(N_permil, probs = 0.25),
            N_q50 = quantile(N_permil, probs = 0.5),
            N_q75 = quantile(N_permil, probs = 0.75),
            ROI_q25 = quantile(ROIAREA, probs = 0.25),
            ROI_q50 = quantile(ROIAREA, probs = 0.5),
            ROI_q75 = quantile(ROIAREA, probs = 0.75),
            n = n()
  )


# PLOT
permil_app_stat_vis <- subset(permil_app_stat, strain=='Marinobacter')

ggplot(permil_app_stat_vis, aes(x=N_q50, y=ROI_q50, shape=ring, colour=treatment)) + 
  geom_point(aes(size=n, fill=treatment, shape=ring), stroke=0.3) + 
  # geom_errorbar(aes(ymax = ROI_q75, ymin = ROI_q25), width=0, linewidth=0.1) +
  # geom_errorbarh(aes(xmax = N_q75, xmin = N_q25), height=0, linewidth=0.1) +
  scale_size(range = c(1, 4)) +
  # scale_x_continuous(transform = 'log') +
  # scale_y_continuous(transform = 'log') +
  # Alcani, Devosi, Marino, none
  scale_color_manual(values=c("#C00000", "#0432FF", "#1a6b3b", "#000000")) +
  scale_fill_manual(values=c("#FFDFE1", "#D2DCFB", "#C4D2CA", "#ffffff")) +
  scale_shape_manual(values = c(16, 21)) +  # Inner, outer
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.2, colour = 'black'),
        axis.ticks = element_line(size = 0.2, colour = 'black'),
        text = element_text(size = 8)
  )

ggsave("figures/fig5b_v5.pdf", width = 1.6, height = 1.6)