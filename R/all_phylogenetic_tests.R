rm(list = ls())
getwd()
setwd("~/GitHub/Phylo_Evol_Timeseries/")

library(phytools)


tree <- read.tree("./data/tree/RAxML_bipartitionsBranchLabels.Task2")

pruned.tree <- drop.tip(tree, c('NC_005042.1.353331-354795', 'Janthinobacterium_sp_KBS0711'))


data <- read.csv(file = './data/pNpS_mean_absolute_differences.csv')

data <- subset(data, taxon!='J')



MAD_0_1 <-  subset(data, treatment_pair == '0_1' )$mean_absolute_difference
names(MAD_0_1) <- subset(data, treatment_pair == '0_1' )$tree_name

MAD_0_2 <-  subset(data, treatment_pair == '0_2' )$mean_absolute_difference
names(MAD_0_2) <- subset(data, treatment_pair == '0_2' )$tree_name

MAD_1_2 <-  subset(data, treatment_pair == '1_2' )$mean_absolute_difference
names(MAD_1_2) <- subset(data, treatment_pair == '1_2' )$tree_name




slopes_se_0_1 <-  subset(data, treatment_pair == '0_1' )$slope_standard_error
names(slopes_se_0_1) <- subset(data, treatment_pair == '0_1' )$tree_name


slopes_0_2 <-  subset(data, treatment_pair == '0_2' )$slope
names(slopes_0_2) <- subset(data, treatment_pair == '0_2' )$tree_name
slopes_se_0_2 <-  subset(data, treatment_pair == '0_2' )$slope_standard_error
names(slopes_se_0_2) <- subset(data, treatment_pair == '0_2' )$tree_name

slopes_1_2 <-  subset(data, treatment_pair == '1_2' )$slope
names(slopes_1_2) <- subset(data, treatment_pair == '1_2' )$tree_name
slopes_se_1_2 <-  subset(data, treatment_pair == '1_2' )$slope_standard_error
names(slopes_se_1_2) <- subset(data, treatment_pair == '1_2' )$tree_name


groups <- c(rep(c('0_1', '0_2', '1_2'), each=5))


MAD_all <- c(MAD_0_1, MAD_0_1)


phylANOVA(pruned.tree, groups, MAD_all)
# 
phyl.pairedttest(pruned.tree, slopes_0_1, x2=slopes_0_2, se1=slopes_0_1, se2=slopes_0_2)

phyl.pairedttest(pruned.tree, slopes_0_1, x2=slopes_1_2, se1=slopes_0_1, se2=slopes_1_2)

phyl.pairedttest(pruned.tree, slopes_0_2, x2=slopes_1_2, se1=slopes_0_2, se2=slopes_1_2)


