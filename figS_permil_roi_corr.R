library("ggplot2")
source('utils.R')

permil_df <- read.csv("data/SIP_permil_v2.csv", check.names = FALSE)
sample_info <- read.csv("data/SIP_sample_info.csv")
sample_info$microplate <- as.character(sample_info$microplate)

permil_append <- append_xnet(permil_df, sample_info)

ggplot(permil_append, aes(x=ROIAREA, y=N_permil)) +
  geom_point(size = 0.6, alpha = 0.4) +
  geom_smooth(method=lm, size = 0.15, alpha = 0.2) +
  facet_wrap(~sample_name, ncol=3, scales = 'free') +
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = 7),
        axis.line = element_line(size = 0.2, colour = 'black'),
        axis.ticks = element_line(colour = 'black', size=0.2),
        # axis.text.y.right = element_blank(),
        # axis.ticks.y.right = element_blank(),
        # axis.line.y.right = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        text = element_text(size = 7)
  )

ggsave("figures/SIP_permil_roi.pdf", width = 4.5, height = 10.5)


ggplot(permil_append, aes(x=C_permil, y=N_permil)) +
  geom_point(size = 0.6, alpha = 0.4) +
  geom_smooth(method=lm, size = 0.15, alpha = 0.2) +
  facet_wrap(~sample_name, ncol=3, scales = 'free') +
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = 7),
        axis.line = element_line(size = 0.2, colour = 'black'),
        axis.ticks = element_line(colour = 'black', size=0.2),
        # axis.text.y.right = element_blank(),
        # axis.ticks.y.right = element_blank(),
        # axis.line.y.right = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        text = element_text(size = 7)
  )

ggsave("figures/SIP_permil_C_N.pdf", width = 4.5, height = 10.5)