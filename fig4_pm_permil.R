## An R-version to analyze incorporation of C and N
library("ggplot2")
library(ggforce)
library("dplyr")
library("ggbreak")
library(FSA)
source('utils.R')


permil_df <- read.csv("data/SIP_permil_v2.csv", check.names = FALSE)
sample_info <- read.csv("data/SIP_sample_info.csv")
sample_info$microplate <- as.character(sample_info$microplate)

permil_append <- append_xnet(permil_df, sample_info)

permil_app_stat <- permil_append %>%
  group_by(treatment, ring) %>%
  summarize(N_q25 = quantile(N_permil, probs = 0.25),
            N_q50 = quantile(N_permil, probs = 0.5),
            N_q75 = quantile(N_permil, probs = 0.75),
            n = n()
  )


# figures for outer strains
permil_append_vis <- subset(permil_append, ring=='outer')
permil_app_stat_vis <- subset(permil_app_stat, ring=='outer')


ggplot() +
  geom_sina(data = permil_append_vis,
            aes(x=treatment, y=N_permil, color=treatment),
            maxwidth = 0.8,
            alpha=0.4,
            size=0.3) +
  geom_errorbar(data=permil_app_stat_vis,
                aes(x=treatment, ymin=N_q25, ymax=N_q75),
                width = 0.15, color='black', size=0.3) +
  geom_errorbar(data=permil_app_stat_vis,
                aes(x=treatment, ymin=N_q50, ymax=N_q50),
                width = 0.3, color='black', size=0.6) +
  # geom_text(data=permil_app_stat_vis,
  #           aes(x=treatment, y=0, label=paste('n = ', n, sep='')),
  #           position=position_dodge(width=0.9), vjust=-0.5, size=2.5) +
  scale_y_continuous(breaks = append(seq(0, 750, 250), seq(1000, 5000, 2000))) +
  scale_y_break(c(750, 751), scales=0.1) +
  # scale_y_continuous(breaks = seq(0, 5000, 250)) +
  scale_color_manual(values=c("#E06666","#5E7BFB","#1a6b3b","#878787")) +
  # Alcani, Devosi, Marino, none
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

ggsave("figures/SIP_N_permil.pdf", width = 1.5, height = 2.5)


# statistical test
kruskal.test(N_permil ~ treatment, data = permil_append_vis)
pairwise.wilcox.test(permil_append_vis$N_permil, permil_append_vis$treatment, p.adjust.method = "BH")

print(
  cat(
    'Alcanivorax:', (df_vis_stat$q50[1]-df_vis_stat$q50[4]) / df_vis_stat$q50[4], 
    '\n',
    'Devosia:', (df_vis_stat$q50[2]-df_vis_stat$q50[4]) / df_vis_stat$q50[4], 
    '\n',
    'Marinobacter:', (df_vis_stat$q50[3]-df_vis_stat$q50[4]) / df_vis_stat$q50[4],
    '\n'
  )
)