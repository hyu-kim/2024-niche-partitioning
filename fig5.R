## Combines count and Cnet data in one plot (day 14)

# import data
count_bact <- read.csv("data/flow2_day14_bact.csv")
count_pt <- read.csv("data/flow2_day14_pt.csv")
cnet_stat <- read.csv("data/SIP_cnet_day14.csv")

# get statistics from count data
