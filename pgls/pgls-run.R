library(ape)
library(geiger)
library(phytools)
library(nlme)
library(scales)
library(gplots)
library(phylolm)

setwd("./pgls")

eggs = read.csv("data/LH_data_reorder.csv",h=T,sep=';')
colnames(eggs)[1] = "species"
colnames(eggs)

eggs$species = gsub(" ", "_", eggs$species)

rownames(eggs) = eggs$species
eggs

### PREPPING THE DATASET ####

eggs$body.log = log(eggs$body.size.avg,10)
eggs$fec.log = log(eggs$fecundity.avg,10)
eggs$eggsize.log = log(eggs$egg.diam.avg,10)

phy = read.tree("data/stern.2017.ml.tree")
plot(phy)
dev.off()

## CLUTCH SIZE TREE 

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


#### PHY SIGNAL ########

### CLUTCH SIZE ###

fec.log = setNames(eggs$fec.log, nm = rownames(eggs))

fec.sig1 = phylosig(pruned, fec.log, test = T, method = "K", nsim = 10^6)
fec.sig2 = phylosig(pruned, fec.log, test = T, method = "lambda", nsim = 10^6)

fec.OU = phylolm(fec.log ~ 1, data = eggs, phy = pruned, model = "OUfixedRoot", 
                   lower.bound = 10e-4, upper.bound = 3)
fec.BM = phylolm(fec.log ~ 1, data = eggs, phy = pruned, model = "BM")

TH_cray = max(branching.times(pruned))

(log(2)/fec.OU$optpar)/TH_cray

fec.BM$aic - fec.OU$aic

## body size correlation ##

fec.OU2 = phylolm(fec.log ~ body.log, data = eggs, phy = pruned, model = "OUfixedRoot", 
                 lower.bound = 10e-4, upper.bound = 5)
fec.BM2 = phylolm(fec.log ~ body.log, data = eggs, phy = pruned, model = "BM")

(log(2)/fec.OU2$optpar)/TH_cray

fec.BM2$aic - fec.OU2$aic

### Egg size ###

diam.log = setNames(diam$eggsize.log, nm = rownames(diam))

diam.sig1 = phylosig(pruned.d, diam.log, test = T, method = "K", nsim = 10^6)
diam.sig2 = phylosig(pruned.d, diam.log, test = T, method = "lambda", nsim = 10^6)

diam.OU = phylolm(diam.log ~ 1, data = diam, phy = pruned.d, model = "OUfixedRoot", 
                 lower.bound = 10e-4, upper.bound = 3)
diam.BM = phylolm(diam.log ~ 1, data = diam, phy = pruned.d, model = "BM")

TH_cray2 = max(branching.times(pruned.d))

(log(2)/diam.OU$optpar)/TH_cray2

diam.BM$aic - diam.OU$aic

## body size correlation ##

diam.OU2 = phylolm(diam.log ~ body.log, data = diam, phy = pruned.d, model = "OUfixedRoot", 
                  lower.bound = 10e-4, upper.bound = 3)
diam.BM2 = phylolm(diam.log ~ body.log, data = diam, phy = pruned.d, model = "BM")


(log(2)/diam.OU2$optpar)/TH_cray2

diam.BM2$aic - diam.OU2$aic



##### OU MODELS FOR CLUTCH SIZE #####

fec.OU.fix = phylolm(fec.log ~ body.log * burrow.lough, data = eggs, phy = pruned, model = "OUfixedRoot", 
                     lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)
fec.OU.rand = phylolm(fec.log ~ body.log * burrow.lough, data = eggs, phy = pruned, model = "OUrandomRoot", 
                     lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)

fec.OU.fix$aic - fec.OU.rand$aic

summary(fec.OU.fix)


plot(fec.log~body.log,data=eggs,cex=1.5,bty='l',las=1,
     pch=21,bg=c('darkorchid1',"orange")[as.numeric(as.factor(eggs$burrow.lough))],
     ylab = "Egg number (log10)", xlab = "Carapace length (log10)")


##### OU MODELS FOR CLUTCH SIZE - TRADITIONAL BURROWING #####

fec.trad.OUF = phylolm(fec.log ~ body.log * traditional.burrow.classification, data = eggs, phy = pruned, model = "OUfixedRoot", 
                     lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)
fec.trad.OUR = phylolm(fec.log ~ body.log * traditional.burrow.classification, data = eggs, phy = pruned, model = "OUrandomRoot", 
                      lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)

fec.trad.OUF$aic - fec.trad.OUR$aic

summary(fec.trad.OUF)


##### OU MODELS FOR CLUTCH SIZE - COMBINED BURROWING #####

fec.comb.OUF = phylolm(fec.log ~ body.log * combined.burrow.classification, data = eggs, phy = pruned, model = "OUfixedRoot", 
                       lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)
fec.comb.OUR = phylolm(fec.log ~ body.log * combined.burrow.classification, data = eggs, phy = pruned, model = "OUrandomRoot", 
                       lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)

fec.comb.OUF$aic - fec.comb.OUR$aic

summary(fec.comb.OUF)
summary(fec.comb.OUR)

##### OU MODELS FOR EGG SIZE #####

diam.OU.fix = phylolm(diam.log ~ body.log * burrow.lough, data = diam, phy = pruned.d, model = "OUfixedRoot", 
                     lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)
diam.OU.rand = phylolm(diam.log ~ body.log * burrow.lough, data = diam, phy = pruned.d, model = "OUrandomRoot", 
                      lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)

diam.OU.fix$aic - diam.OU.rand$aic

summary(diam.OU.fix)


##### OU MODELS FOR EGG SIZE - TRADITIONAL BURROWING #####

diam.trad.OUF = phylolm(diam.log ~ body.log * traditional.burrow.classification, data = diam, phy = pruned.d, model = "OUfixedRoot", 
                      lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)
diam.trad.OUR = phylolm(diam.log ~ body.log * traditional.burrow.classification, data = diam, phy = pruned.d, model = "OUrandomRoot", 
                       lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)

diam.trad.OUF$aic - diam.trad.OUR$aic

summary(diam.trad.OUF)
summary(diam.trad.OUR)


##### OU MODELS FOR EGG SIZE - TRADITIONAL BURROWING #####

diam.comb.OUF = phylolm(diam.log ~ body.log * combined.burrow.classification, data = diam, phy = pruned.d, model = "OUfixedRoot", 
                        lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)
diam.comb.OUR = phylolm(diam.log ~ body.log * combined.burrow.classification, data = diam, phy = pruned.d, model = "OUrandomRoot", 
                        lower.bound = 10e-4, upper.bound = 10, boot = 10^2, measurement_error = T)

diam.comb.OUF$aic - diam.comb.OUR$aic

summary(diam.comb.OUF)
summary(diam.comb.OUR)

############ ANCESTRAL STATE BURROWING ##############

### MAXIMUM LIKELIHOOD ###
burrow = setNames(eggs$burrow.lough,rownames(eggs))
burrow.trad = setNames(eggs$traditional.burrow.classification,rownames(eggs))
burrow.comb = setNames(eggs$combined.burrow.classification,rownames(eggs))

bur.ER = fitDiscrete(pruned, burrow, model = "ER", control = list(niter = 1000), ncores = 4)

bur.ARD = fitDiscrete(pruned, burrow, model = "ARD", control = list(niter = 1000), ncores = 4)

bur.ER$opt$aicc
bur.ARD$opt$aicc

bur.trad.ER = fitDiscrete(pruned, burrow.trad, model = "ER", control = list(niter = 1000), ncores = 4)

bur.trad.ARD = fitDiscrete(pruned, burrow.trad, model = "ARD", control = list(niter = 1000), ncores = 4)

bur.trad.ER$opt$aicc
bur.trad.ARD$opt$aicc

saveRDS(bur.trad.ER,"bur-trad-ER.RDS")
saveRDS(bur.trad.ARD,"bur-trad-ARD.RDS")

bur.comb.ER = fitDiscrete(pruned, burrow.comb, model = "ER", control = list(niter = 1000), ncores = 4)

bur.comb.ARD = fitDiscrete(pruned, burrow.comb, model = "ARD", control = list(niter = 1000), ncores = 4)

bur.comb.ER$opt$aicc
bur.comb.ARD$opt$aicc

saveRDS(bur.comb.ER,"bur-comb-ER.RDS")
saveRDS(bur.comb.ARD,"bur-comb-ARD.RDS")

### They have similar AICcs. So, ER model wins because it is simpler.

#bur.sim = make.simmap(pruned, burrow, model = "ER", nsim = 1000, Q = "empirical")
#bur.sim2 = make.simmap(pruned, burrow, model = "ER", nsim = 1000, Q = "mcmc")

#saveRDS(bur.sim, file = 'data/burrow-simmap1.RDS')
#saveRDS(bur.sim2, file = 'data/burrow-simmap2.RDS')

bur.sim = readRDS(file = 'data/burrow-simmap1.RDS')
bur.sim2 = readRDS(file = 'data/burrow-simmap2.RDS')

summary(bur.sim)
summary(bur.sim2)

plot(density(bur.sim))

png("figures/simmap-morph.png", res = 600, units = 'mm',
    w = 260, h = 220)
plot(density(bur.sim2))
dev.off()

mean(logLik(bur.sim))
mean(logLik(bur.sim2))

#bur.sim is better, has less degrees of freedom and slightly lower log likelihood


pdf("figures/simmap.pdf", h = 20, w = 20)
plot(summary(bur.sim), type = 'fan', ftype = 'i')
nodelabels(piecol=cols2,cex=0.22)

dev.off()

pdf("figures/simmap2.pdf", h = 20, w = 20)
plot(summary(bur.sim2), type = 'fan', ftype = 'i')
dev.off()


#### ANC_RECON TRADITIONAL ####

#bur.trad = make.simmap(pruned, burrow.trad, model = "ER", nsim = 1000, Q = "empirical")
#bur.trad2 = make.simmap(pruned, burrow.trad, model = "ER", nsim = 1000, Q = "mcmc")

#saveRDS(bur.trad, file = 'data/burrow-trad-simmap.RDS')
#saveRDS(bur.trad2, file = 'data/burrow-trad-simmap2.RDS')

bur.trad = readRDS(file = 'data/burrow-trad-simmap.RDS')
bur.trad2 = readRDS(file = 'data/burrow-trad-simmap2.RDS')

summary(bur.trad)
summary(bur.trad2)

png("figures/simmap-trad.png", res = 600, units = 'mm',
    w = 260, h = 220)
par(mfrow=c(3,2))
plot(density(bur.trad), transition = 'primary_burrower->secondary_burrower')
legend("topleft","a)",bty='n')
plot(density(bur.trad), transition = 'secondary_burrower->tertiary_burrower')
legend("topleft","b)",bty='n')
plot(density(bur.trad), transition = 'primary_burrower->tertiary_burrower')
legend("topleft","c)",bty='n')

plot(density(bur.trad), transition = 'secondary_burrower->primary_burrower')
legend("topleft","d)",bty='n')
plot(density(bur.trad), transition = 'tertiary_burrower->primary_burrower')
legend("topleft","e)",bty='n')
dev.off()

plot(density(bur.trad2))

mean(logLik(bur.trad))
mean(logLik(bur.trad2))

#bur.trad is better, has less degrees of freedom and slightly lower log likelihood

pdf("figures/simmap-trad.pdf", h = 20, w = 20)
plot(summary(bur.trad), type = 'fan', ftype = 'i')
nodelabels(piecol=cols2,cex=0.22)

dev.off()

pdf("figures/simmap2-trad.pdf", h = 20, w = 20)
plot(summary(bur.trad2), type = 'fan', ftype = 'i')
dev.off()


#### ANC_RECON COMBINED ####

#bur.comb = make.simmap(pruned, burrow.comb, model = "ER", nsim = 1000, Q = "empirical")
#bur.comb2 = make.simmap(pruned, burrow.comb, model = "ER", nsim = 1000, Q = "mcmc")

#saveRDS(bur.comb, file = 'data/burrow-comb-simmap.RDS')
#saveRDS(bur.comb2, file = 'data/burrow-comb-simmap2.RDS')

bur.comb = readRDS(file = 'data/burrow-comb-simmap.RDS')
bur.comb2 = readRDS(file = 'data/burrow-comb-simmap2.RDS')

summary(bur.comb)
summary(bur.comb2)

png("figures/simmap-comb.png", res = 600, units = 'mm',
    w = 260, h = 220)
par(mfrow=c(1,1))
plot(density(bur.comb))
dev.off()

plot(density(bur.comb2))

mean(logLik(bur.comb))
mean(logLik(bur.comb2))


pdf("figures/simmap-comb.pdf", h = 20, w = 20)
plot(summary(bur.comb), type = 'fan', ftype = 'i')
nodelabels(piecol=cols2,cex=0.22)

dev.off()

pdf("figures/simmap2-comb.pdf", h = 20, w = 20)
plot(summary(bur.comb2), type = 'fan', ftype = 'i')
dev.off()
