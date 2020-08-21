rm(list = ls())
getwd()
setwd("~/GitHub/Phylo_Evol_Timeseries/")

set.seed(123456)

library('ape')


# Load ML tree
ml.tree <- read.tree("data/tree/RAxML_bipartitionsBranchLabels.Task2")
# Define the outgroup
outgroup <- match("NC_005042.1.353331-354795", ml.tree$tip.label)
# Create a rooted tree {ape} using the outgroup
ml.rooted <- root(ml.tree, outgroup, resolve.root = TRUE)
# Keep rooted but drop outgroup branch
ml.rooted <- drop.tip(ml.rooted, c("NC_005042.1.353331-354795"))
is.ultrametric(ml.rooted)
ml.rooted.um  <- chronos(ml.rooted)
is.ultrametric(ml.rooted.um)

plot(ml.rooted.um)

prop.part(ml.rooted.um)
phy <- reorder(ml.rooted.um, "postorder")
phy$tip.label <- c("B", "D", "F", "P", "J", "C")
C<-vcv.phylo(ml.rooted.um)

# test data matrix
col1.1 <- rnorm(2, mean = 0, sd = 1)
col1.2 <- rnorm(2, mean = 3, sd = 0.3)
col2.1 <- rnorm(2, mean = 0, sd = 1)
col2.2 <- rnorm(2, mean = 3, sd = 0.3)
Y <- data.frame(col1 = c(col1.1,col1.2), col2= c(col2.1,col2.2))
rownames(Y) <- c("0C1_100", "0D1_100", "0C1_200", "0D1_200")

# get new covariance matrix
C.new <- as.data.frame(matrix(ncol=nrow(Y), nrow=nrow(Y)))
rownames(C.new) <- rownames(Y)
colnames(C.new) <- rownames(Y)

for (row_i in rownames(C.new))
{
  for (row_j in rownames(C.new))
  {
    taxon_i <- row_i[2]
    taxon_j <- row_j[2]
    print(taxon_i)
    print(taxon_j)
  }
}


# function to compute phylogenetic VCV using joint Pagel's lambda
# written by Liam Revell 2011
phyl.vcv<-function(X,phy,lambda){
  C<-vcv.phylo(phy)
  C<-lambda.transform(lambda,C)
  invC<-solve(C)
  a<-matrix(colSums(invC%*%X)/sum(invC),ncol(X),1)
  A<-matrix(rep(a,nrow(X)),nrow(X),ncol(X),byrow=T)
  V<-t(X-A)%*%invC%*%(X-A)/(nrow(C)-1)
  return(list(C=C,R=V,alpha=a))
}





# Brownian motion model
temp<-phyl.vcv(Y,C,1.0)
V<-temp$R
a<-t(temp$alpha)
C<-temp$C






