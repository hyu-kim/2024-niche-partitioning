## An R-version to analyze flow cytometry reads
## Compute a fold-increase, then correct by algae number
library("Rmisc")
library("ggplot2")
library(ggforce)

df <- read.csv(file='data/flow2_bact.csv')
df <- df[df$Time!=0,]

df_alg <- read.csv(file='data/flow2_pt.csv') # algae

# replace treatment values
df$Treatment[df$Treatment==1] <- 'Devosia'
df$Treatment[df$Treatment==2] <- 'none'
df$Treatment[df$Treatment==3] <- 'Alcanivorax'
df$Treatment[df$Treatment==4] <- 'Marinobacter'
df_alg$Treatment[df_alg$Treatment==1] <- 'Devosia'
df_alg$Treatment[df_alg$Treatment==2] <- 'none'
df_alg$Treatment[df_alg$Treatment==3] <- 'Alcanivorax'
df_alg$Treatment[df_alg$Treatment==4] <- 'Marinobacter'



## Compute fold increase, day 5 -- 14
df_fold <- data.frame(matrix(nrow=0, ncol=7))
colnames(df_fold) <- c('Strain', 'Treatment', 'Microplate', 'Ring', 'Direction',
                       'fold', 'fold_adj')

conds <- c('Treatment', 'Microplate', 'Ring', 'Direction')
df_cond <- df[,conds]

for(i in 1:nrow(unique(df_cond))){
  ind <- (df$Treatment==df_cond[i,1]) & (df$Microplate==df_cond[i,2]) & 
    (df$Ring==df_cond[i,3]) & (df$Direction==df_cond[i,4])
  df_ext <- df[ind,]
  
  ind2 <- (df_alg$Treatment==df_cond[i,1]) & (df_alg$Microplate==df_cond[i,2])
  df_alg_ext <- df_alg[ind2,]
  if (length(df_ext[,1])>1){
    val <- df_ext[2,8] / df_ext[1,8]
    val2 <- df_alg_ext[2,5] / df_alg_ext[1,5]
    df_app <- data.frame(Strain=df_ext[1,2], 
                         Treatment=df_ext[1,4], 
                         Microplate=df_ext[1,5], 
                         Ring=df_ext[1,6], 
                         Direction=df_ext[1,7], 
                         fold=val,
                         fold_adj=val/val2)
    df_fold[nrow(df_fold)+1, ] = df_app
  }
}

df_fold$Strain <- factor(df_fold$Strain, 
                         levels=c('Devosia EAB7WZ', 'Alcanivorax EA2', 
                                  'None', 'Marinobacter 3-2'))
df_fold$Treatment <- factor(df_fold$Treatment, 
                            levels=c('Alcanivorax', 'Devosia',  
                                     'Marinobacter', 'none')
                            )
df_fold$fold_log <- log(df_fold$fold)  # natural log
df_fold$fold_adj_log <- log(df_fold$fold_adj)



## subset abundance on day 14 only
df14 <- df[df$Time==14,]
df_alg14 <- df_alg[df_alg$Time==14,]
df14$Abd_per_alga <- matrix(0, nrow=nrow(df14))
conds <- c('Treatment', 'Microplate')
df14_cond <- df14[,conds]
df14_cond_unique <- unique(df14_cond)

# export day14
write.csv(df14, "data/flow2_day14_bact.csv", row.names=FALSE)
write.csv(df_alg14, "data/flow2_day14_pt.csv", row.names=FALSE)

for(i in 1:nrow(df14_cond_unique)){
  ind <- (df14$Treatment==df14_cond_unique[i,1]) & 
    (df14$Microplate==df14_cond_unique[i,2])
  
  ind2 <- (df_alg14$Treatment==df14_cond_unique[i,1]) & 
    (df_alg14$Microplate==df14_cond_unique[i,2])
  
  df14[ind, 9] <- df14[ind, 8] / df_alg14[ind2, 5]
}
df14$Strain <- factor(df14$Strain, 
                      levels=c('Alcanivorax EA2', 'Devosia EAB7WZ',
                               'Marinobacter 3-2', 'None'))
df14$Treatment <- factor(df14$Treatment, levels=c(3,1,4,2))
df14_se <- summarySE(df14, measurevar = 'Abd_per_alga', groupvars = c('Treatment', 'Ring'))



## subset day 5
df5 <- df[df$Time==5,]
df5$Strain <- factor(df5$Strain,
                     levels=c('Alcanivorax EA2', 'Devosia EAB7WZ',
                              'Marinobacter 3-2', 'None'))
df5$Treatment <- factor(df5$Treatment, levels=c(3,1,4,2))