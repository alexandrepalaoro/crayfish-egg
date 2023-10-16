#### This code was used to read the output from the mcmc and steppingstone runs.
#### We ran both mcmcs and steppingstones in a cluster separately, and analyzed them here.
#### The parsing_code.R and egg_parsing_code.R are very similar codes because they do similar things.
#### THe difference will rely in the models with the best fit and the figures made.

require(bayou)

## Setting workdirectory
setwd("./egg_size")

## Making sure is okay
getwd()

m#### DATA DRIVEN MODELS ####

#### 11 ####
### Reading the results from the mcmc model in which slope and intercept were global (i.e., they did not vary)

mcmc.11_r1 = readRDS("result_mcmc/mcmc11-RUN1.rds")
mcmc.11_r2 = readRDS("result_mcmc/mcmc11-RUN2.rds")
mcmc.11_r3 = readRDS("result_mcmc/mcmc11-RUN3.rds")
mcmc.11_r4 = readRDS("result_mcmc/mcmc11-RUN4.rds")

chain.11_r1 <- set.burnin(mcmc.11_r1$load(), 0.3)
chain.11_r2 <- set.burnin(mcmc.11_r2$load(), 0.3)
chain.11_r3 <- set.burnin(mcmc.11_r3$load(), 0.3)
chain.11_r4 <- set.burnin(mcmc.11_r4$load(), 0.3)

### This won't work because there aren't any shifts. I'm just making sure it reflects what is expected.
shiftsums11_r1 <- shiftSummaries(chain.11_r1, mcmc.11_r1, pp.cutoff=0.7)
shiftsums11_r2 <- shiftSummaries(chain.11_r2, mcmc.11_r2, pp.cutoff=0.7)
shiftsums11_r3 <- shiftSummaries(chain.11_r3, mcmc.11_r3, pp.cutoff=0.7)
shiftsums11_r4 <- shiftSummaries(chain.11_r4, mcmc.11_r4, pp.cutoff=0.7)

### Combining the chains that converged. 
### The Gelman plots that show model convergence are located in the code used to run the model in the cluster.

r11_mcmc.chain.list <- list(chain.11_r1, chain.11_r2, chain.11_r3, chain.11_r4)
r11.mcmc <- combine.chains(chain.list = r11_mcmc.chain.list)

saveRDS(r11.mcmc, "result_mcmc/combined_11.RDS")

#### N1 ####
### Reading the results from the mcmc model in which did not vary slope but intercept could vary.

mcmc.N1_r1 = readRDS("result_mcmc/mcmcN1-RUN1.rds")
mcmc.N1_r2 = readRDS("result_mcmc/mcmcN1-RUN2.rds")
mcmc.N1_r3 = readRDS("result_mcmc/mcmcN1-RUN3.rds")
mcmc.N1_r4 = readRDS("result_mcmc/mcmcN1-RUN4.rds")

chain.N1_r1 <- set.burnin(mcmc.N1_r1$load(), 0.3)
chain.N1_r2 <- set.burnin(mcmc.N1_r2$load(), 0.3)
chain.N1_r3 <- set.burnin(mcmc.N1_r3$load(), 0.3)
chain.N1_r4 <- set.burnin(mcmc.N1_r4$load(), 0.3)

summary(chain.N1_r1)

### Checking the shifts, just to have an idea if the chains give the same results.
shiftsumsN1_r1 <- shiftSummaries(chain.N1_r1, mcmc.N1_r1, pp.cutoff=0.7)
shiftsumsN1_r2 <- shiftSummaries(chain.N1_r2, mcmc.N1_r2, pp.cutoff=0.7)
shiftsumsN1_r3 <- shiftSummaries(chain.N1_r3, mcmc.N1_r3, pp.cutoff=0.7)
shiftsumsN1_r4 <- shiftSummaries(chain.N1_r4, mcmc.N1_r4, pp.cutoff=0.7)

# Checking the descendents
shiftsumsN1_r1$descendents
shiftsumsN1_r2$descendents
shiftsumsN1_r3$descendents
shiftsumsN1_r4$descendents

# Ploting them
plotShiftSummaries(shiftsumsN1_r1, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsN1_r2, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsN1_r3, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsN1_r4, single.plot=TRUE, label.pts = FALSE)

### Combining the chains that converged. 
### The Gelman plots that show model convergence are located in the code used to run the model in the cluster.
N1_mcmc.chain.list <- list(chain.N1_r2, chain.N1_r3)
N1.mcmc <- combine.chains(chain.list = N1_mcmc.chain.list)

saveRDS(N1.mcmc, "result_mcmc/combined_N1.RDS")

#### NN ####
### Reading the results from the mcmc model in which both slope but intercept could vary.

mcmc.NN_r1 = readRDS("result_mcmc/mcmcNN-RUN1.rds")
mcmc.NN_r2 = readRDS("result_mcmc/mcmcNN-RUN2.rds")
mcmc.NN_r3 = readRDS("result_mcmc/mcmcNN-RUN3.rds")
mcmc.NN_r4 = readRDS("result_mcmc/mcmcNN-RUN4.rds")

chain.NN_r1 <- set.burnin(mcmc.NN_r1$load(), 0.3)
chain.NN_r2 <- set.burnin(mcmc.NN_r2$load(), 0.3)
chain.NN_r3 <- set.burnin(mcmc.NN_r3$load(), 0.3)
chain.NN_r4 <- set.burnin(mcmc.NN_r4$load(), 0.3)

### Checking the shifts, just to have an idea if the chains give the same results.
shiftsumsNN_r1 <- shiftSummaries(chain.NN_r1, mcmc.NN_r1, pp.cutoff=0.7)
shiftsumsNN_r2 <- shiftSummaries(chain.NN_r2, mcmc.NN_r2, pp.cutoff=0.7)
shiftsumsNN_r3 <- shiftSummaries(chain.NN_r3, mcmc.NN_r3, pp.cutoff=0.7)
shiftsumsNN_r4 <- shiftSummaries(chain.NN_r4, mcmc.NN_r4, pp.cutoff=0.7)

# Checking the descendents
shiftsumsNN_r1$descendents
shiftsumsNN_r2$descendents
shiftsumsNN_r3$descendents
shiftsumsNN_r4$descendents

# Ploting them
plotShiftSummaries(shiftsumsNN_r1, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsNN_r2, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsNN_r3, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsNN_r4, single.plot=TRUE, label.pts = FALSE)

### Combining the chains that converged. 
### The Gelman plots that show model convergence are located in the code used to run the model in the cluster.
NN_mcmc.chain.list <- list(chain.NN_r1, chain.NN_r2, chain.NN_r3, chain.NN_r4)
NN.mcmc <- combine.chains(chain.list = NN_mcmc.chain.list)

saveRDS(NN.mcmc, "result_mcmc/combined_NN.RDS")

#### MORPHOLOGICAL CLASSIFICATION #####
### In this model, we fixed the locations where shifts could occur and  allowed both the intercept and slope to vary.

mcmc.F_r1 = readRDS("result_mcmc/mcmcF-RUN1.rds")
mcmc.F_r2 = readRDS("result_mcmc/mcmcF-RUN2.rds")
mcmc.F_r3 = readRDS("result_mcmc/mcmcF-RUN3.rds")
mcmc.F_r4 = readRDS("result_mcmc/mcmcF-RUN4.rds")

chain.F_r1 <- set.burnin(mcmc.F_r1$load(), 0.3)
chain.F_r2 <- set.burnin(mcmc.F_r2$load(), 0.3)
chain.F_r3 <- set.burnin(mcmc.F_r3$load(), 0.3)
chain.F_r4 <- set.burnin(mcmc.F_r4$load(), 0.3)

### Checking the shifts, just to have an idea if the chains give the same results.
shiftsumsF_r1 <- shiftSummaries(chain.F2_r1, mcmc.F_r1, pp.cutoff=0.7)
shiftsumsF_r2 <- shiftSummaries(chain.F2_r2, mcmc.F_r2, pp.cutoff=0.7)
shiftsumsF_r3 <- shiftSummaries(chain.F2_r3, mcmc.F_r3, pp.cutoff=0.7)
shiftsumsF_r4 <- shiftSummaries(chain.F2_r4, mcmc.F_r4, pp.cutoff=0.7)

# Checking the descendents
shiftsumsF_r1$descendents
shiftsumsF_r2$descendents
shiftsumsF_r3$descendents
shiftsumsF_r4$descendents

# Ploting them
plotShiftSummaries(shiftsumsF_r1, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsF_r2, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsF_r3, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsF_r4, single.plot=TRUE, label.pts = FALSE)

### Combining the chains that converged. 
### The Gelman plots that show model convergence are located in the code used to run the model in the cluster.
F_mcmc.chain.list <- list(chain.F_r1, chain.F_r3, chain.F_r4)
F.mcmc <- combine.chains(chain.list = F_mcmc.chain.list)

saveRDS(F.mcmc, "result_mcmc/combined_F.RDS")

#### COMBINED CLASSIFICATION #####
### In this model, we fixed the locations where shifts could occur and  allowed both the intercept and slope to vary.

mcmc.C_r1 = readRDS("result_mcmc/mcmcC-RUN1.rds")
mcmc.C_r2 = readRDS("result_mcmc/mcmcC-RUN2.rds")
mcmc.C_r3 = readRDS("result_mcmc/mcmcC-RUN3.rds")
mcmc.C_r4 = readRDS("result_mcmc/mcmcC-RUN4.rds")

chain.C_r1 <- set.burnin(mcmc.C_r1$load(), 0.3)
chain.C_r2 <- set.burnin(mcmc.C_r2$load(), 0.3)
chain.C_r3 <- set.burnin(mcmc.C_r3$load(), 0.3)
chain.C_r4 <- set.burnin(mcmc.C_r4$load(), 0.3)

### Checking the shifts, just to have an idea if the chains give the same results.
shiftsumsC_r1 <- shiftSummaries(chain.C_r1, mcmc.C_r1, pp.cutoff=0.7)
shiftsumsC_r2 <- shiftSummaries(chain.C_r2, mcmc.C_r2, pp.cutoff=0.7)
shiftsumsC_r3 <- shiftSummaries(chain.C_r3, mcmc.C_r3, pp.cutoff=0.7)
shiftsumsC_r4 <- shiftSummaries(chain.C_r4, mcmc.C_r4, pp.cutoff=0.7)

# Checking the descendents
shiftsumsC_r1$descendents
shiftsumsC_r2$descendents
shiftsumsC_r3$descendents
shiftsumsC_r4$descendents

# Ploting them
plotShiftSummaries(shiftsumsC_r1, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsC_r2, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsC_r3, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsC_r4, single.plot=TRUE, label.pts = FALSE)

### Combining the chains that converged. 
### The Gelman plots that show model convergence are located in the code used to run the model in the cluster.
C_mcmc.chain.list <- list(chain.C_r1, chain.C_r2, chain.C_r3, chain.C_r4)
C.mcmc <- combine.chains(chain.list = C_mcmc.chain.list)

saveRDS(C.mcmc, "result_mcmc/combined_C.RDS")

#### TRADITIONAL CLASSIFICATION #####
### In this model, we fixed the locations where shifts could occur and  allowed both the intercept and slope to vary.

mcmc.T_r1 = readRDS("result_mcmc/mcmcT-RUN1.rds")
mcmc.T_r2 = readRDS("result_mcmc/mcmcT-RUN2.rds")
mcmc.T_r3 = readRDS("result_mcmc/mcmcT-RUN3.rds")
mcmc.T_r4 = readRDS("result_mcmc/mcmcT-RUN4.rds")

chain.T_r1 <- set.burnin(mcmc.T_r1$load(), 0.3)
chain.T_r2 <- set.burnin(mcmc.T_r2$load(), 0.3)
chain.T_r3 <- set.burnin(mcmc.T_r3$load(), 0.3)
chain.T_r4 <- set.burnin(mcmc.T_r4$load(), 0.3)

### Checking the shifts, just to have an idea if the chains give the same results.
shiftsumsT_r1 <- shiftSummaries(chain.T_r1, mcmc.T_r1, pp.cutoff=0.7)
shiftsumsT_r2 <- shiftSummaries(chain.T_r2, mcmc.T_r2, pp.cutoff=0.7)
shiftsumsT_r3 <- shiftSummaries(chain.T_r3, mcmc.T_r3, pp.cutoff=0.7)
shiftsumsT_r4 <- shiftSummaries(chain.T_r4, mcmc.T_r4, pp.cutoff=0.7)

# Checking the descendents
shiftsumsT_r1$descendents
shiftsumsT_r2$descendents
shiftsumsT_r3$descendents
shiftsumsT_r4$descendents

# Ploting them
plotShiftSummaries(shiftsumsT_r1, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsT_r2, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsT_r3, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsT_r4, single.plot=TRUE, label.pts = FALSE)

### Combining the chains that converged. 
### The Gelman plots that show model convergence are located in the code used to run the model in the cluster.
T_mcmc.chain.list <- list(chain.T_r1, chain.T_r3, chain.T_r4)
T.mcmc <- combine.chains(chain.list = T_mcmc.chain.list)

saveRDS(T.mcmc, "result_mcmc/combined_T.RDS")



####### Checking shifts and the differences between the shifts in the models ########
## To do this, we need to load the combined chains, but we also needs the model used to run the mcmc. NOT the individual chain, but the 
## commands and arguments used to build the mcmc chain. That is why I'm also loading one of the mcmcs.

mcmc.11_r1 = readRDS("result_mcmc/mcmc11-RUN1.rds")
chain.11 = readRDS("result_mcmc/combined_11.RDS") 

mcmc.N1_r1 = readRDS("result_mcmc/mcmcN1-RUN1.rds")
chain.N1 = readRDS("result_mcmc/combined_N1.RDS") 

mcmc.NN_r1 = readRDS("result_mcmc/mcmcNN-RUN1.rds")
chain.NN = readRDS("result_mcmc/combined_NN.RDS") 

mcmc.F_r1 = readRDS("result_mcmc/mcmcF-RUN1.rds")
chain.F = readRDS("result_mcmc/combined_F.RDS") 

mcmc.C_r1 = readRDS("result_mcmc/mcmcC-RUN1.rds")
chain.C = readRDS("result_mcmc/combined_C.RDS") 

mcmc.T_r1 = readRDS("result_mcmc/mcmcT-RUN1.rds")
chain.T = readRDS("result_mcmc/combined_T.RDS") 

## Loading the inputed variables for the models, like the tree, the data and body length. We will need that in some cases.
pruned = mcmc.NN_r1$tree
dat = mcmc.NN_r1$dat
eggs = mcmc.NN_r1$pred

### Now, let's build the shift summaries for each model using two posterior probabilit cutoffs 
shiftsumsN1 = shiftSummaries(chain.N1, mcmc.N1_r1, pp.cutoff=0.7)
shiftsumsN1.5 = shiftSummaries(chain.N1, mcmc.N1_r1, pp.cutoff=0.5)

shiftsumsNN = shiftSummaries(chain.NN, mcmc.NN_r1, pp.cutoff=0.7)
shiftsumsNN.5 = shiftSummaries(chain.NN, mcmc.NN_r1, pp.cutoff=0.5)

shiftsumsF = shiftSummaries(chain.F, mcmc.F_r1, pp.cutoff=0.7)
shiftsumsF.5 = shiftSummaries(chain.F, mcmc.F_r1, pp.cutoff=0.5)

shiftsumsC = shiftSummaries(chain.C, mcmc.C_r1, pp.cutoff=0.7)
shiftsumsC.5 = shiftSummaries(chain.C, mcmc.C_r1, pp.cutoff=0.5)

shiftsumsT = shiftSummaries(chain.T, mcmc.T_r1, pp.cutoff=0.7)
shiftsumsT.5 = shiftSummaries(chain.T, mcmc.T_r1, pp.cutoff=0.5)

#ploting shifts
plotShiftSummaries(shiftsumsN1, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsN1.5, lwd=2, single.plot=TRUE, label.pts=FALSE)

plotShiftSummaries(shiftsumsNN, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsNN.5, lwd=2, single.plot=TRUE, label.pts=FALSE)

plotShiftSummaries(shiftsumsF, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsF.5, lwd=2, single.plot=TRUE, label.pts=FALSE)

plotShiftSummaries(shiftsumsC, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsC.5, lwd=2, single.plot=TRUE, label.pts=FALSE)

plotShiftSummaries(shiftsumsT, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsT.5, lwd=2, single.plot=TRUE, label.pts=FALSE)


plotSimmap.mcmc(chain.N1, burnin=0.3, pp.cutoff=0.7)
plotSimmap.mcmc(chain.N1, burnin=0.3, pp.cutoff=0.5)

plotSimmap.mcmc(chain.NN, burnin=0.3, pp.cutoff=0.7)
plotSimmap.mcmc(chain.NN, burnin=0.3, pp.cutoff=0.5)

plotSimmap.mcmc(chain.F, burnin=0.3, pp.cutoff = 0.7)
plotSimmap.mcmc(chain.F, burnin=0.3, pp.cutoff = 0.5)

plotSimmap.mcmc(chain.C, burnin=0.3, pp.cutoff = 0.7)
plotSimmap.mcmc(chain.C, burnin=0.3, pp.cutoff = 0.5)

plotSimmap.mcmc(chain.T, burnin=0.3, pp.cutoff = 0.7)
plotSimmap.mcmc(chain.T, burnin=0.3, pp.cutoff = 0.5)


phenogram.density(pruned,dat,burnin=0.3,chain.11,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.11,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.N1,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.N1,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.NN,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.NN,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.F,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.F,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.C,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.C,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.T,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.T,pp.cutoff = 0.5)


#### Now we can check how each parameter shifts across the tree.
## For the egg size models, the best model had shifts in both intercept and sloped. Thus, we're plotting all of them to check them out.

plotBranchHeatMap(pruned, chain.N1, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))

plotBranchHeatMap(pruned, chain.NN, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))
plotBranchHeatMap(pruned, chain.NN, variable="beta_logBody", burnin=0.3, pal=rainbow, legend_ticks=seq(0.2, 1.5, 0.1))

plotBranchHeatMap(pruned, chain.F, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))
plotBranchHeatMap(pruned, chain.F, variable="beta_logBody", burnin=0.3, pal=rainbow, legend_ticks=seq(0.2, 1.5, 0.1))

plotBranchHeatMap(pruned, chain.C, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))
plotBranchHeatMap(pruned, chain.C, variable="beta_logBody", burnin=0.3, pal=rainbow, legend_ticks=seq(0.2, 1.5, 0.1))

plotBranchHeatMap(pruned, chain.T, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))
plotBranchHeatMap(pruned, chain.T, variable="beta_logBody", burnin=0.3, pal=rainbow, legend_ticks=seq(0.2, 1.5, 0.1))


#### STEPSTONES ####

### Here, we will load the stepstones andcalculate the bayes factor to know which model is better.

ss.11 = readRDS("result_mcmc/ss_11.rds")
ss.N1 = readRDS("result_mcmc/ss_N1.rds")
ss.NN = readRDS("result_mcmc/ss_NN.rds")
ss.F = readRDS("result_mcmc/ss_F.rds")
ss.C = readRDS("result_mcmc/ss_C.rds")
ss.T = readRDS("result_mcmc/ss_T.rds")


mlnL <- c("11"=ss.11$lnr, "N1"=ss.N1$lnr, "NN"=ss.NN$lnr, "M"=ss.F$lnr, "C"=ss.C$lnr, "T" = ss.T$lnr)
mlnL


summ_BF <- as.data.frame(c(2*(mlnL[1]-mlnL[3]), 2*(mlnL[2]-mlnL[3]), 2*(mlnL[4]-mlnL[3]),
                           2*(mlnL[5]-mlnL[3]), 2*(mlnL[6]-mlnL[3])))
rownames(summ_BF) <- c("NN vs 11", "NN vs N1", "NN vs M", "NN vs C", "NN vs T")
colnames(summ_BF) <- c("Bayes_Factor")
summ_BF


##### FIGURES FOR THE PAPER ######

#This first one is directly related to the bayes factor.

png("figures/eggsize-Bayesfactor.png", res = 600, unit = 'mm', h = 160, w = 160)
plot(summ_BF$Bayes_Factor, xlab="Models", ylab="Bayes Factor", pch=19,
     xaxt = 'n', cex = 2.2, bty = 'l', las = 1, ylim = c(0,100))
axis(side = 1, at = c(1:5), labels = rownames(summ_BF))
abline(h=10, lwd = 2, lty = 2, col = 'red')
dev.off()


## Now, the figures with the shift sums.
# Given NN was the best model, all figures will be based on it.

mcmc.NN = readRDS("result_mcmc/combined_NN.RDS")
mcmc.NN_r1 = readRDS("result_mcmc/mcmcNN-RUN1.rds")

pruned = mcmc.NN_r1$tree
dat = mcmc.NN_r1$dat
body = mcmc.NN_r1$pred

plot(dat ~ body)

#plot(logFec ~ logBody, data = eggs)

shiftsums.NN <- shiftSummaries(mcmc.NN, mcmc.NN_r1, pp.cutoff=0.5)
plotShiftSummaries(shiftsums.NN, single.plot=TRUE, label.pts = FALSE)

shiftsums.NN$descendents


pdf("figures/simmap.eggsize-0.5pcutoff.pdf",h=18,w=12)
plotSimmap.mcmc(mcmc.NN, burnin=0.3, edge.type = "regimes", pp.cutoff=0.5, no.margin=TRUE,
                pal=colorRampPalette(colors=c("#E0B27F","#6B2DE0")))
dev.off()

shiftsums.NN$descendents
shiftsums.NN$regressions
shiftsums.NN$cladesummaries[[1]]$densities

estimates = shiftsums.NN$regressions


eggs = read.csv("data/LH_data.csv",h=T,sep=';')
colnames(eggs)[1] = "species"
eggs$species = gsub(" ", "_", eggs$species)

eggs$logBody = log(eggs$body.size.avg)
eggs$logFec = log(eggs$fecundity.avg)

camb = eggs[eggs$genus=="Cambarellus",]

camb$logEgg = log(camb$egg.diam.avg,10)
camb$logBody = log(camb$body.size.avg,10)

row.names(dat) = eggs$species
dat = as.matrix(dat)[,1]


logBody = as.matrix(data.frame(eggs$logBody)) 
rownames(logBody) = eggs$species

row.names(dat) = eggs$species
dat = as.matrix(dat)[,1]

png("figures/carapace-egg-pp0.5.png",res=600,units='mm',h=150,w=200)

tiff("figures/carapace-egg-pp0.5.tiff",res=600,units='mm',h=150,w=200,
    compression = 'lzw')
plot(dat ~ body, cex = 2, las = 1, bty = 'l', 
     ylab = "Egg size (log10)", xlab = "Body length (log10)",
     pch = 21, bg = "#E0B27F")


points(x = camb$logBody, y = camb$logEgg, cex = 2, pch = 21, bg = "#6B2DE0")

curve(estimates[1,1]+(estimates[1,2]*x), add =T, lwd = 3, col = "#E0B27F" )
curve(estimates[2,1]+(estimates[2,2]*x), add =T, lwd = 3, col = "#6B2DE0", from = 1, to =1.3)
dev.off()