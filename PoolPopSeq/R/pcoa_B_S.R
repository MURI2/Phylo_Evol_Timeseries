rm(list = ls())
getwd()
setwd("~/GitHub/Task2/PoolPopSeq/")


package.list <- c('vegan', 'dplyr', 'tidyr', 'BiodiversityR', 'indicspecies', 'viridis')
for (package in package.list) {
  if (!require(package, character.only=TRUE, quietly=TRUE)) {
    install.packages(package)
    library(package, character.only=TRUE)
  }
}


######
# run PCoA for gene_by_pop containing G scores for all time points
gene_by_pop <- c("data/gene_by_sample/B_S/sample_by_gene_Gscore.txt")
df.all <- read.table(paste(gene_by_pop, collapse = ''), sep = "\t", header = TRUE, row.names = 1)
df.all.no0 <- df.all[rowSums(df.all[,-1]) != 0,]
df.all.db <- vegdist(df.all.no0, method = "bray", upper = TRUE, diag = TRUE)
df.all.pcoa <- cmdscale(df.all.db, eig = TRUE, k = 14) 
explainvar1 <- round(df.all.pcoa$eig[1] / sum(df.all.pcoa$eig), 3) * 100
explainvar2 <- round(df.all.pcoa$eig[2] / sum(df.all.pcoa$eig), 3) * 100
explainvar3 <- round(df.all.pcoa$eig[3] / sum(df.all.pcoa$eig), 3) * 100
sum.eig <- sum(explainvar1, explainvar2, explainvar3)
png(filename = paste(c("figs/pcoa/pcoa_B_S.png"), collapse = ''),
    width = 1200, height = 1200, res = 96*2)
# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)
# Initiate Plot
plot(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2], xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
cols <- c()
treats <- c()
for (x in rownames(df.all.pcoa$points)){
  if (grepl("frequency_L0", x)){
    treats <- c(treats, "1")
    cols <- c(cols, "#87CEEB")
  } else if ( grepl("frequency_L1", x)) {
    treats <- c(treats, "10")
    cols <- c(cols, "#FFA500")
  } else if (grepl("frequency_L2", x)) {
    treats <- c(treats, "100")
    cols <- c(cols, "#FF6347")
  }
}

times <- c()
for (x in rownames(df.all.pcoa$points)){
  if (grepl("D100", x)){
    times <- c(times, "100")
  } else if ( grepl("D200", x)) {
    times <- c(times, "200")
  } else if (grepl("D300", x)) {
    times <- c(times, "300")
  }
}

#plot
points(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2],
       pch = 19, cex = 3, bg = "gray", col = cols)
ordiellipse(df.all.pcoa, treats, conf = 0.95)
text(df.all.pcoa$points[ ,1], df.all.pcoa$points[ ,2], labels=times, cex= 0.7)
dev.off()




# Define Order of Sites
order <- rev(attr(df.all.db.100, "Labels"))  

# Plot Heatmap
levelplot(as.matrix(df.all.db.100)[, order], aspect = "iso", col.regions = inferno, 
          xlab = "Doubs Site", ylab = "Doubs Site", scales = list(cex = 0.5), 
          main = "Bray-Curtis Distance", type="lower")




######
# run PCoA for gene_by_pop containing G scores for just day 100
df.all.no0.100 <- df.all.no0[grep('D100', rownames(df.all.no0)), ]
df.all.db.100 <- vegdist(df.all.no0.100, method = "bray", upper = TRUE, diag = TRUE)
df.all.pcoa.100 <- cmdscale(df.all.db.100, eig = TRUE, k = 2) 
explainvar1 <- round(df.all.pcoa.100$eig[1] / sum(df.all.pcoa.100$eig), 3) * 100
explainvar2 <- round(df.all.pcoa.100$eig[2] / sum(df.all.pcoa.100$eig), 3) * 100
#explainvar3 <- round(df.all.pcoa.100$eig[3] / sum(df.all.pcoa.100$eig), 3) * 100
png(filename = paste(c("figs/pcoa/pcoa_D100_B_S.png"), collapse = ''),
    width = 1200, height = 1200, res = 96*2)
# Define Plot Parameters
#par(mar = c(5, 5, 1, 2) + 0.1)
B_S.D100 <- c()
for (x in rownames(df.all.pcoa.100$points)){
  if (grepl("B", x)){
    B_S.D100 <- c(B_S.D100, 19)
  } else if ( grepl("S", x)) {
    B_S.D100 <- c(B_S.D100, 1)
  } 
}
par(mar = c(6.5, 6, 1.5, 2.5) + 0.1)
# Initiate Plot
print(df.all.pcoa.100$points[ ,1])
print(df.all.pcoa.100$points[ ,2])

plot(df.all.pcoa.100$points[ ,1], df.all.pcoa.100$points[ ,2],  xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
cols.D100 <- c()
treats.D100 <- c()
for (x in rownames(df.all.pcoa.100$points)){
  if (grepl("frequency_L0", x)){
    treats.D100 <- c(treats.D100, "1")
    cols.D100 <- c(cols.D100, "#87CEEB")
  } else if ( grepl("frequency_L1", x)) {
    treats.D100 <- c(treats.D100, "10")
    cols.D100 <- c(cols.D100, "#FFA500")
  } else if (grepl("frequency_L2", x)) {
    treats.D100 <- c(treats.D100, "100")
    cols.D100 <- c(cols.D100, "#FF6347")
  }
}
test <- cbind.data.frame(B_S.D100, treats.D100)
test$merge <- paste(test$B_S.D100, test$treats.D100, sep="_")
#plot
points(df.all.pcoa.100$points[ ,1], df.all.pcoa.100$points[ ,2],
       pch = B_S.D100, cex = 2.5, bg = "gray", col = cols.D100, lwd  = 3)
ellipse.cols <- c("#87CEEB", "#FFA500", "#FF6347", "#87CEEB", "#FFA500", "#FF6347")
ordiellipse(df.all.pcoa.100, test$merge, conf = 0.95, col = ellipse.cols, lwd=2)
legend(x=-0.8,y = -1.01, xpd = TRUE, c("1-day","10-day","100-day"), 
       col=c('#87CEEB','#FFA500','#FF6347'), ncol=3, bty ="n", 
       pch = c(16,16,16), cex = 1.6, pt.cex = 3.0, text.font = 20)
dev.off()



########
# Just the ellipses 
png(filename = paste(c("figs/pcoa/pcoa_D100_B_S_ellipses.png"), collapse = ''),
    width = 1200, height = 1200, res = 96*2)
# Define Plot Parameters
par(mar = c(6.5, 6, 1.5, 2.5) + 0.1)
# Initiate Plot
plot(df.all.pcoa.100$points[ ,1], df.all.pcoa.100$points[ ,2],  xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
#plot
ordiellipse(df.all.pcoa.100, test$merge, conf = 0.95, col = ellipse.cols, lwd=2)
#text(df.all.pcoa.100$points[ ,1], df.all.pcoa.100$points[ ,2], labels=B_S.D100, cex= 0.7)
legend(x=-0.8,y = -1.01, xpd = TRUE, c("1-day","10-day","100-day"), 
       col=c('#87CEEB','#FFA500','#FF6347'), ncol=3, bty ="n", 
       pch = c(16,16,16), cex = 1.6, pt.cex = 3.0, text.font = 20)
dev.off()



#######
# plot for 1, 10, and 100 days lines
df.all.no0.100.1day.points <- df.all.pcoa.100$points[grep('L0', rownames(df.all.pcoa.100$points)), ]
png(filename = paste(c("figs/pcoa/pcoa_D100_B_S_1day.png"), collapse = ''),
    width = 1200, height = 1200, res = 96*2)
# Define Plot Parameters
B_S.D100.1day <- c()
for (x in rownames(df.all.no0.100.1day.points)){
  if (grepl("B", x)){
    B_S.D100.1day <- c(B_S.D100.1day, 19)
  } else if ( grepl("S", x)) {
    B_S.D100.1day <- c(B_S.D100.1day, 1)
  } 
}
par(mar = c(6.5, 6, 1.5, 2.5) + 0.1)
# Initiate Plot
plot(df.all.no0.100.1day.points[ ,1], df.all.no0.100.1day.points[ ,2], xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
cols.D100.1day <- c()
treats.D100.1day <- c()
for (x in rownames(df.all.no0.100.1day.points)){
  if (grepl("frequency_L0", x)){
    treats.D100.1day <- c(treats.D100.1day, "1")
    cols.D100.1day <- c(cols.D100.1day, "#87CEEB")
  } else if ( grepl("frequency_L1", x)) {
    treats.D100.1day <- c(treats.D100.1day, "10")
    cols.D100.1day <- c(cols.D100.1day, "#FFA500")
  } else if (grepl("frequency_L2", x)) {
    treats.D100.1day <- c(treats.D100.1day, "100")
    cols.D100.1day <- c(cols.D100.1day, "#FF6347")
  }
}
#plot
points(df.all.no0.100.1day.points[ ,1], df.all.no0.100.1day.points[ ,2],
       pch = B_S.D100.1day, cex = 2.5, bg = "gray", col = cols.D100.1day, lwd  = 3)
legend(x=-0.8,y = -1.01, xpd = TRUE, c("1-day","10-day","100-day"), 
       col=c('#87CEEB','#FFA500','#FF6347'), ncol=3, bty ="n", 
       pch = c(16,16,16), cex = 1.6, pt.cex = 3.0, text.font = 20)
dev.off()
#######
####10 day
######
df.all.no0.100.10day.points <- df.all.pcoa.100$points[grep('L1', rownames(df.all.pcoa.100$points)), ]
png(filename = paste(c("figs/pcoa/pcoa_D100_B_S_10day.png"), collapse = ''),
    width = 1200, height = 1200, res = 96*2)
# Define Plot Parameters
B_S.D100.10day <- c()
for (x in rownames(df.all.no0.100.10day.points)){
  if (grepl("B", x)){
    B_S.D100.10day <- c(B_S.D100.10day, 19)
  } else if ( grepl("S", x)) {
    B_S.D100.10day <- c(B_S.D100.10day, 1)
  } 
}
par(mar = c(6.5, 6, 1.5, 2.5) + 0.1)
# Initiate Plot
plot(df.all.no0.100.10day.points[ ,1], df.all.no0.100.10day.points[ ,2], xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
cols.D100.10day <- c()
treats.D100.10day <- c()
for (x in rownames(df.all.no0.100.10day.points)){
  if (grepl("frequency_L0", x)){
    treats.D100.10day <- c(treats.D100.10day, "1")
    cols.D100.10day <- c(cols.D100.10day, "#87CEEB")
  } else if ( grepl("frequency_L1", x)) {
    treats.D100.10day <- c(treats.D100.10day, "10")
    cols.D100.10day <- c(cols.D100.10day, "#FFA500")
  } else if (grepl("frequency_L2", x)) {
    treats.D100.10day <- c(treats.D100.10day, "100")
    cols.D100.10day <- c(cols.D100.10day, "#FF6347")
  }
}
#plot
points(df.all.no0.100.10day.points[ ,1], df.all.no0.100.10day.points[ ,2],
       pch = B_S.D100.10day, cex = 2.5, bg = "gray", col = cols.D100.10day, lwd  = 3)
legend(x=-0.8,y = -1.01, xpd = TRUE, c("1-day","10-day","100-day"), 
       col=c('#87CEEB','#FFA500','#FF6347'), ncol=3, bty ="n", 
       pch = c(16,16,16), cex = 1.6, pt.cex = 3.0, text.font = 20)
dev.off()
#######
####100 day
######
df.all.no0.100.100day.points <- df.all.pcoa.100$points[grep('L2', rownames(df.all.pcoa.100$points)), ]
png(filename = paste(c("figs/pcoa/pcoa_D100_B_S_100day.png"), collapse = ''),
    width = 1200, height = 1200, res = 96*2)
# Define Plot Parameters
B_S.D100.100day <- c()
for (x in rownames(df.all.no0.100.100day.points)){
  if (grepl("B", x)){
    B_S.D100.100day <- c(B_S.D100.100day, 19)
  } else if ( grepl("S", x)) {
    B_S.D100.100day <- c(B_S.D100.100day, 1)
  } 
}
par(mar = c(6.5, 6, 1.5, 2.5) + 0.1)
# Initiate Plot
plot(df.all.no0.100.100day.points[ ,1], df.all.no0.100.100day.points[ ,2],  xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
cols.D100.100day <- c()
treats.D100.100day <- c()
for (x in rownames(df.all.no0.100.100day.points)){
  if (grepl("frequency_L0", x)){
    treats.D100.100day <- c(treats.D100.100day, "1")
    cols.D100.100day <- c(cols.D100.100day, "#87CEEB")
  } else if ( grepl("frequency_L1", x)) {
    treats.D100.100day <- c(treats.D100.100day, "10")
    cols.D100.100day <- c(cols.D100.100day, "#FFA500")
  } else if (grepl("frequency_L2", x)) {
    treats.D100.100day <- c(treats.D100.100day, "100")
    cols.D100.100day <- c(cols.D100.100day, "#FF6347")
  }
}
#plot
points(df.all.no0.100.100day.points[ ,1], df.all.no0.100.100day.points[ ,2],
       pch = B_S.D100.100day, cex = 2.5, bg = "gray", col = cols.D100.100day, lwd  = 3)
legend(x=-0.8,y = -1.01, xpd = TRUE, c("1-day","10-day","100-day"), 
       col=c('#87CEEB','#FFA500','#FF6347'), ncol=3, bty ="n", 
       pch = c(16,16,16), cex = 1.6, pt.cex = 3.0, text.font = 20)
#legend(x=-0.7,y = -1.01, xpd = TRUE, c("1-day","10-day","100-day"), 
#       col=c('#87CEEB','#FFA500','#FF6347'), ncol=3, bty ="n", 
#       pch = c(16,16,16), cex = 1.6, pt.cex = 3.0, text.font = 20)
dev.off()







#########
#PERMANOVA
#########
B_S.D100.category <- c()
for (x in rownames(df.all.no0.100)){
  if (grepl("B", x)){
    B_S.D100.category <- c(B_S.D100.category, 'B')
  } else if ( grepl("S", x)) {
    B_S.D100.category <- c(B_S.D100.category, 'S')
  } 
}

treats.D100.perm <- c()
for (x in rownames(df.all.no0.100)){
  if (grepl("frequency_L0", x)){
    treats.D100.perm <- c(treats.D100.perm, 0)
  } else if ( grepl("frequency_L1", x)) {
    treats.D100.perm <- c(treats.D100.perm, 1)
  } else if (grepl("frequency_L2", x)) {
    treats.D100.perm <- c(treats.D100.perm, 2)
  }
}


perm <- adonis(df.all.no0.100 ~ B_S.D100.category * treats.D100.perm, method = "bray", permutations = 999)

interact.beta <- perm$coef.sites[4,]
dim(t(perm$coef.sites) *  perm$model.matrix)
t(perm$coef.sites) *  perm$model.matrix
perm$coef.sites * t(perm$model.matrix)
write.table(t(perm$coef.sites) *  perm$model.matrix, "data/betas.txt", sep="\t")

# plot coefficients 









# Add Points & Labels
treats.D100.indval <- c()
for (x in rownames(df.all.no0.100)){
  if (grepl("frequency_L0", x)){
    treats.D100.indval <- c(treats.D100.indval, "1")
  } else if ( grepl("frequency_L1", x)) {
    treats.D100.indval <- c(treats.D100.indval, "10")
  } else if (grepl("frequency_L2", x)) {
    treats.D100.indval <- c(treats.D100.indval, "100")
  }
}

B_S.D100.indval <- c()
for (x in rownames(df.all.no0.100)){
  if (grepl("B", x)){
    B_S.D100.indval <- c(B_S.D100.indval, "W")
  } else if ( grepl("S", x)) {
    B_S.D100.indval <- c(B_S.D100.indval, "S")
  } 
}

test.indval <- cbind.data.frame(B_S.D100.indval, treats.D100.indval)
test.indval$merge <- paste(test.indval$B_S.D100, test.indval$treats.D100, sep="_")

indval <- multipatt(df.all.no0.100, cluster = test.indval$merge, func = "IndVal.g", control = how(nperm=999)) 
summary(indval)



if(strain == 'S'){
  days <- c('D100')
} else {
  days <- c('D100', 'D200', 'D300')
}
for (day in days){
  to.keep <- c()
  df.all.no0.day <- df.all.no0[grep(day, rownames(df.all.no0)), ]
  treats.day <- c()
  for (x in rownames(df.all.no0.day)){
    if (grepl("frequency_L0", x)){
      treats.day <- c(treats.day, "1")
    } else if ( grepl("frequency_L1", x)) {
      treats.day <- c(treats.day, "10")
    } else if (grepl("frequency_L2", x)) {
      treats.day <- c(treats.day, "100")
    }
  }
  #print(day)
  # Run PERMANOVA with adonis function
  #print(adonis(df.all.no0.day ~ treats.day, method = "bray", permutations = 999))
  # Indicator Value
  #indval <- multipatt(df.all.no0.day, cluster = treats.day, func = "IndVal.g", control = how(nperm=1000)) 
  #summary(indval)
}
# calculate distance between vector time-points 
transfer.times <- c('0', '1', '2')
reps <- c('1', '2', '3', '4', '5')
days <- c(100, 200)
fileName <- paste("data/euclidean_distance/", strain, ".txt", sep = "")
sink(fileName)
header <- paste("Line","Time1", "Time2", "Distance", sep = "\t")
cat(header)
cat("\n")
for (transfer.time in transfer.times){
  for (rep in reps){
    for (day in days){
      v1 <- paste("frequency_L", transfer.time, strain, rep, '_D', toString(day), sep = "")
      v2 <- paste("frequency_L", transfer.time, strain, rep, '_D', toString(day+100), sep = "")
      v1.exist <- any(row.names(df.all.pcoa$points) == v1)
      v2.exist <- any(row.names(df.all.pcoa$points) == v2)
      if (v1.exist && v2.exist){
        dist.v1.v2 <- dist(rbind(df.all.pcoa$points[v1,], df.all.pcoa$points[v2,]))[1]
        line.out <- paste(paste("L", transfer.time, strain, rep, sep="") , toString(day), toString(day+100), toString(dist.v1.v2), sep = "\t")
        cat(line.out)
        cat("\n")
      }
    }
  }
}
sink()
