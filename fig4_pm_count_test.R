# Continued from "fig4_pm_count.R"

## Perform T-test on fold increase
df_fold2 <- df_fold[df_fold$Ring==2,]
df_fold1 <- df_fold[df_fold$Ring==1,]
df_fold_se <- summarySE(df_fold, measurevar = "fold", groupvars = c("Strain","Treatment","Ring"))
df_fold_se2 <- summarySE(df_fold, measurevar = "fold2_log", groupvars = c("Strain","Treatment","Ring"))
# no correction (bacteria only)
t.test(df_fold2[df_fold2$Treatment==4,8], df_fold2[df_fold2$Treatment==1,8])
t.test(df_fold2[df_fold2$Treatment==4,8], df_fold2[df_fold2$Treatment==2,8])
t.test(df_fold2[df_fold2$Treatment==4,8], df_fold2[df_fold2$Treatment==3,8])
# corrected by algal number
t.test(df_fold2[df_fold2$Treatment==2,9], df_fold2[df_fold2$Treatment==1,9])
t.test(df_fold2[df_fold2$Treatment==1,9], df_fold2[df_fold2$Treatment==3,9])
t.test(df_fold2[df_fold2$Treatment==4,9], df_fold2[df_fold2$Treatment==3,9])



## Perform t-test for each timepoint (4: Marino, 3: Alcani, 1: Devosi, 2: Empty)
t.test(df5[(df5$Treatment==4)&(df5$Ring==2),8], df5[(df5$Treatment==1)&(df5$Ring==2),8])
t.test(df5[(df5$Treatment==4)&(df5$Ring==2),8], df5[(df5$Treatment==2)&(df5$Ring==2),8])
t.test(df5[(df5$Treatment==4)&(df5$Ring==2),8], df5[(df5$Treatment==3)&(df5$Ring==2),8])
t.test(df5[(df5$Treatment==3)&(df5$Ring==2),8], df5[(df5$Treatment==1)&(df5$Ring==2),8])
t.test(df5[(df5$Treatment==3)&(df5$Ring==2),8], df5[(df5$Treatment==2)&(df5$Ring==2),8])
t.test(df5[(df5$Treatment==1)&(df5$Ring==2),8], df5[(df5$Treatment==2)&(df5$Ring==2),8])

t.test(df14[(df14$Treatment==4)&(df14$Ring==2),8], df14[(df14$Treatment==1)&(df14$Ring==2),8])
t.test(df14[(df14$Treatment==4)&(df14$Ring==2),8], df14[(df14$Treatment==2)&(df14$Ring==2),8])
t.test(df14[(df14$Treatment==4)&(df14$Ring==2),8], df14[(df14$Treatment==3)&(df14$Ring==2),8])
t.test(df14[(df14$Treatment==3)&(df14$Ring==2),8], df14[(df14$Treatment==1)&(df14$Ring==2),8])
t.test(df14[(df14$Treatment==3)&(df14$Ring==2),8], df14[(df14$Treatment==2)&(df14$Ring==2),8])
t.test(df14[(df14$Treatment==1)&(df14$Ring==2),8], df14[(df14$Treatment==2)&(df14$Ring==2),8])

t.test(df[(df$Treatment==3)&(df$Ring==2),8], df[(df$Treatment==2)&(df$Ring==2),8])