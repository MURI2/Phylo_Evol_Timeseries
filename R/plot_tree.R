rm(list = ls())
getwd()
setwd("~/GitHub/Phylo_Evol_Timeseries/")

library('seqinr')
library('ggtree')
library('treeio')


ml.tree <- ggtree::read.raxml("data/tree/RAxML_bipartitionsBranchLabels.Task2")

to.keep <- ml.tree@phylo$tip.label
to.keep <- to.keep[to.keep!= 'NC_005042.1.353331-354795']
pruned.tree <- treeio::drop.tip(ml.tree, ml.tree@phylo$tip.label[-match(to.keep, ml.tree@phylo$tip.label)])


tree.plot <- ggtree(pruned.tree, branch.length='none') +
  #geom_text(aes(label=bootstrap)) +
  geom_tiplab(size = 3) + 
  geom_label2(aes(label = bootstrap), size=2.5, label.size = 0.3) +
  geom_treescale(x = 0, y = -0.5) + xlim(0, 10)


ggsave(file="figs/tree.png", tree.plot, units='in', dpi=600)
