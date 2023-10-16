library(ape)
library(geiger)
library(phytools)
library(nlme)
library(scales)
library(gplots)
library(phylolm)
library(plotrix)

setwd(./pgls")

eggs <- read.csv("LH_data_reorder.csv",h=T,sep=';')
colnames(eggs)[1] = "species"
colnames(eggs)

eggs$species = gsub(" ", "_", eggs$species)

rownames(eggs) = eggs$species
eggs

### PREPPING THE DATASET ####

eggs$body.log = log(eggs$body.size.avg,10)
eggs$fec.log = log(eggs$fecundity.avg,10)
eggs$eggsize.log = log(eggs$egg.diam.avg,10)

phy = read.tree("stern.2017.ml.tree")
plot(phy)
dev.off()

## EGG NUMBER TREE 

pruned = drop.tip(phy,phy$tip.label[-match(eggs$species,phy$tip.label)])

plot(pruned,no.margin=T)
dev.off()

name.check(pruned,eggs)

## EGG SIZE TREE

diam <- subset(eggs,egg.diam.avg > 0.01)
colnames(diam)
rownames(diam)

pruned.d = drop.tip(phy,phy$tip.label[-match(diam$species,phy$tip.label)])



plot(pruned.d,no.margin=T)
dev.off()

name.check(pruned.d,diam)

############ ANCESTRAL STATE BURROWING ##############
### MAXIMUM LIKELIHOOD ###
##Morphological Classification (listed here as burrow.lough)
burrow = setNames(eggs$burrow.lough,rownames(eggs))

bur.ER = fitDiscrete(pruned, burrow, model = "ER", control = list(niter = 1000), ncores = 4)

bur.ARD = fitDiscrete(pruned, burrow, model = "ARD", control = list(niter = 1000), ncores = 4)

bur.ER$opt$aicc
bur.ARD$opt$aicc

### They have similar AICcs. So, ER model wins because it is simpler.

bur.sim = make.simmap(pruned, burrow, model = "ER", nsim = 1000, Q = "empirical")

eggs = read.csv("LH_data_reorder.csv",h=T,sep=';')


############################ HERE IS ALL YOU NEED TO DO TO MAKE THE PLOT #####################################
t <- max(nodeHeights(pruned))
t #262MYA

test <- summary(bur.sim, plot = FALSE)
test

colorbars2 <-setNames(c("#bc672c", "#79c5e7"), c("Burrower", "Non-burrower"))
ss2 <-getStates(pruned, "tips")
str(ss2)
barcols2 <- setNames(sapply(ss2, function(x, y) y[which(names(y) == x)], y = colorbars2), names(ss2))

colors2 <- c("#bc672c", "#79c5e7")
cols2 <- setNames(colors2[1:length(unique(x))], sort(unique(x)))


#png("Burrow_phylo_with_tip_label.png",res=300,w=200,h=160,units='mm')

plotTree(pruned, type = "fan", fsize = 0.5, ftype = "i", offset = 3, part = 0.98)

tiplabels(pie = to.matrix(burrow[pruned$tip.label], unique(sort(burrow))), piecol = cols2, cex = 0.20)
nodelabels(pie = test$ace, piecol = cols2, cex = 0.3)
add.simmap.legend(colors=cols2,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(pruned)),fsize=0.8)

tick.spacing <- 20
min.tick <- 0
scale <- axis(1, pos = -5, at = seq(t, min.tick, by = -tick.spacing), cex.axis = 0.5, labels = FALSE)

for(i in 1:length(scale)){
  a1 <- 0
  a2 <- 2*pi
  draw.arc(0, 0, radius = scale[i], a1, a2, lwd = 1,
           col = make.transparent("grey", 0.15))
}
text(scale, rep(-23, length(scale)), t - scale, cex = 0.5)
text(mean(scale), -38, "time (mya)", cex = 0.5)

#dev.off()
#############################################################################################################3
#png("Burrow_phylo_without_tip_label.png",res=300,w=200,h=160,units='mm')

plotTree(pruned, type = "fan", fsize = 0.5, ftype = "off", offset = 3, part = 0.98)

tiplabels(pie = to.matrix(burrow[pruned$tip.label], unique(sort(burrow))), piecol = cols2, cex = 0.30)
nodelabels(pie = test$ace, piecol = cols2, cex = 0.45)
add.simmap.legend(colors=cols2,prompt=FALSE,x=0.8*par()$usr[1],
                  y=-max(nodeHeights(pruned)),fsize=0.8)

tick.spacing <- 20
min.tick <- 0
scale <- axis(1, pos = -5, at = seq(t, min.tick, by = -tick.spacing), cex.axis = 0.5, labels = FALSE)

for(i in 1:length(scale)){
  a1 <- 0
  a2 <- 2*pi
  draw.arc(0, 0, radius = scale[i], a1, a2, lwd = 1,
           col = make.transparent("grey", 0.15))
}
text(scale, rep(-23, length(scale)), t - scale, cex = 0.5)
text(mean(scale), -38, "time (mya)", cex = 0.5)

#dev.off()