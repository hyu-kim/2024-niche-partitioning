## Combines count, Cnet, biovolume in one plot (day 14)
library("dplyr")
library(ggplot2)
library(scales)

# import data
count_bact <- read.csv("data/flow2_day14_bact.csv")
count_pt <- read.csv("data/flow2_day14_pt.csv")
cnet <- read.csv("data/SIP_cnet.csv")
cnet_info <- read.csv("data/SIP_sample_info.csv")
mass <- read.table("data/mass.csv", header = TRUE, sep = ',')  # [fg]

count_bact$Ring[count_bact$Ring==1] <- 'inner'
count_bact$Ring[count_bact$Ring==2] <- 'outer'

# Here we assume that bacterial cell number and Cnet are conditionally 
# independent given algal cell number. For each microplate, the amount of 
# bacterial incorporation of algal carbon is estimated by multiplying Cnet median
# with total bacterial number per ring. It is then divided by algal cell number
# to disentangle the algal effect. Normalized values are used for statistical
# analysis, e.g., average or sd.

# summarize count data
count_bact_summ <- count_bact %>%
  group_by(Treatment, Ring, Microplate) %>%
  summarize(abd_mean = mean(Abundance), 
            abd_sd = sd(Abundance),
            abd_med = median(Abundance),
            count_total = 6 * median(Abundance)
            ) # merge over 6 wells per ring

# summarize cnet
cnet$microplate <- NA
cnet$ring <- NA
cnet$treatment <- NA
cnet$strain <- NA
for (s in unique(cnet_info$sample_name)){
  m <- cnet_info$microplate[cnet_info$sample_name==s]
  r <- cnet_info$ring[cnet_info$sample_name==s]
  t <- cnet_info$treatment[cnet_info$sample_name==s]
  st <- cnet_info$strain[cnet_info$sample_name==s]
  
  cnet$microplate[cnet$sample_name==s] <- m
  cnet$ring[cnet$sample_name==s] <- r
  cnet$treatment[cnet$sample_name==s] <- t
  cnet$strain[cnet$sample_name==s] <- st
}

cnet_summ <- cnet %>%
  group_by(treatment, microplate, ring, strain) %>%
  summarize(cnet_q50 = quantile(Cnet, probs = 0.5), 
            cnet_q25 = quantile(Cnet, probs = 0.25),
            cnet_q75 = quantile(Cnet, probs = 0.75)
            )

# merge count, cnet, mass, algal count in to df
df <- select(cnet_summ, treatment, microplate, ring, strain, cnet_q50)
df$count_total <- NA
df$mass_c <- NA
df$count_pt <- NA
df <- df[df$strain!='none',]

for (row in 1:nrow(df)){
  t <- df$treatment[row]
  m <- df$microplate[row]
  r <- df$ring[row]
  st <- df$strain[row]
  df[row, 'count_total'] <- count_bact_summ$count_total[
    count_bact_summ$Treatment == t & 
      count_bact_summ$Ring == r &
      count_bact_summ$Microplate == m
    ]
  df[row, 'mass_c'] <- mass$mass[mass$Genus==st]
  df[row, 'count_pt'] <- count_pt$Abundance[
    count_pt$Treatment==t & 
      count_pt$Microplate==m
    ]
}

df$c_incorp <- df$cnet_q50 * df$count_total * df$mass_c
df$c_incorp_pt <- df$c_incorp / (df$count_pt * mass[4, 'mass'])

# remove Alcani, microplate 2 from df, confirming with XM
# df <- df[!(df$treatment=='Alcanivorax' & df$microplate==2),]

# plot
pdf("figures/fig5b.pdf", width = 2.5, height = 2)

ggplot(df, aes(x=microplate, 
               y=c_incorp_pt, 
               fill=treatment, 
               alpha=factor(ring, levels=c("outer","inner"))
               )
       ) + 
  geom_bar(stat = 'identity', position = 'stack', color='black', size=0.2) +
  facet_grid(~factor(treatment, 
                     levels=c('Devosia', 'Marinobacter', 'Alcanivorax', 'none')
                     )) +
  # geom_text(aes(label=round(c_incorp_pt, digits = 4)), vjust=0.5, alpha=1) +
  # ylim(0, 0.08) +
  # ylab("Incorporated algal carbon (um3)") +
  scale_fill_manual(values=c("#E06666", "#4061f4", "#AB7942", "#5b5b5b")) +  # Alcani, Devosi, Marino, none
  scale_alpha_discrete(range=c(0.3, 1)) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank(),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(size = 0.2),
        strip.text.x = element_blank()
        # axis.ticks = element_blank()
        )

dev.off()