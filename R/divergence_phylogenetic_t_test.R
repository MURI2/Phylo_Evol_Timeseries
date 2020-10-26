rm(list = ls())
getwd()
setwd("~/GitHub/Phylo_Evol_Timeseries/")

library(phytools)


tree <- read.tree("./data/tree/RAxML_bipartitionsBranchLabels.Task2")

pruned.tree <- drop.tip(tree, c('NC_005042.1.353331-354795', 'Janthinobacterium_sp_KBS0711'))


data <- read.csv(file = './data/divergence_slopes.csv')

data <- subset(data, taxon!='J')

slopes_0_1 <-  subset(data, treatment_pair == '0_1' )$slope
names(slopes_0_1) <- subset(data, treatment_pair == '0_1' )$tree_name
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




# 
phyl.pairedttest(pruned.tree, slopes_0_1, x2=slopes_0_2, se1=slopes_0_1, se2=slopes_0_2)

phyl.pairedttest(pruned.tree, slopes_0_1, x2=slopes_1_2, se1=slopes_0_1, se2=slopes_1_2)

phyl.pairedttest(pruned.tree, slopes_0_2, x2=slopes_1_2, se1=slopes_0_2, se2=slopes_1_2)


