#### This code was used to read the output from the mcmc and steppingstone runs.
#### We ran both mcmcs and steppingstones in a cluster separately, and analyzed them here.
#### The parsing_code.R and egg_parsing_code.R are very similar codes because they do similar things.
#### THe difference will rely in the models with the best fit and the figures made.

require(bayou)

## Setting workdirectory
setwd("./cluster")

## Making sure is okay
getwd()


#### DATA DRIVEN MODELS ####

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
N1_mcmc.chain.list <- list(chain.N1_r2, chain.N1_r3, chain.N1_r4)
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
NN_mcmc.chain.list <- list(chain.NN_r1, chain.NN_r2, chain.NN_r3)
NN.mcmc <- combine.chains(chain.list = NN_mcmc.chain.list)

saveRDS(NN.mcmc, "result_mcmc/combined_NN.RDS")

#### MORPHOLOGICAL CLASSIFICATION #####
### In this model, the fixed the locations where shifts could occur and only allowed the intercept to vary. Slope was global

mcmc.F2_r1 = readRDS("result_mcmc/mcmcF2-RUN1.rds")
mcmc.F2_r2 = readRDS("result_mcmc/mcmcF2-RUN2.rds")
mcmc.F2_r3 = readRDS("result_mcmc/mcmcF2-RUN3.rds")
mcmc.F2_r4 = readRDS("result_mcmc/mcmcF2-RUN4.rds")

chain.F2_r1 <- set.burnin(mcmc.F2_r1$load(), 0.3)
chain.F2_r2 <- set.burnin(mcmc.F2_r2$load(), 0.3)
chain.F2_r3 <- set.burnin(mcmc.F2_r3$load(), 0.3)
chain.F2_r4 <- set.burnin(mcmc.F2_r4$load(), 0.3)

### Checking the shifts, just to have an idea if the chains give the same results.
shiftsumsF2_r1 <- shiftSummaries(chain.F2_r1, mcmc.F2_r1, pp.cutoff=0.7)
shiftsumsF2_r2 <- shiftSummaries(chain.F2_r2, mcmc.F2_r2, pp.cutoff=0.7)
shiftsumsF2_r3 <- shiftSummaries(chain.F2_r3, mcmc.F2_r3, pp.cutoff=0.7)
shiftsumsF2_r4 <- shiftSummaries(chain.F2_r4, mcmc.F2_r4, pp.cutoff=0.7)

# Checking the descendents
shiftsumsF2_r1$descendents
shiftsumsF2_r2$descendents
shiftsumsF2_r3$descendents
shiftsumsF2_r4$descendents

# Ploting them
plotShiftSummaries(shiftsumsF2_r1, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsF2_r2, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsF2_r3, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsF2_r4, single.plot=TRUE, label.pts = FALSE)

### Combining the chains that converged. 
### The Gelman plots that show model convergence are located in the code used to run the model in the cluster.
F2_mcmc.chain.list <- list(chain.F2_r1, chain.F2_r2, chain.F2_r3, chain.F2_r4)
F2.mcmc <- combine.chains(chain.list = F2_mcmc.chain.list)

saveRDS(F2.mcmc, "result_mcmc/combined_F2.RDS")

#### COMBINED CLASSIFICATION #####
### In this model, the fixed the locations where shifts could occur and only allowed the intercept to vary. Slope was global

mcmc.C2_r1 = readRDS("result_mcmc/mcmcC2-RUN1.rds")
mcmc.C2_r2 = readRDS("result_mcmc/mcmcC2-RUN2.rds")
mcmc.C2_r3 = readRDS("result_mcmc/mcmcC2-RUN3.rds")
mcmc.C2_r4 = readRDS("result_mcmc/mcmcC2-RUN4.rds")

chain.C2_r1 <- set.burnin(mcmc.C2_r1$load(), 0.3)
chain.C2_r2 <- set.burnin(mcmc.C2_r2$load(), 0.3)
chain.C2_r3 <- set.burnin(mcmc.C2_r3$load(), 0.3)
chain.C2_r4 <- set.burnin(mcmc.C2_r4$load(), 0.3)

### Checking the shifts, just to have an idea if the chains give the same results.
shiftsumsC2_r1 <- shiftSummaries(chain.C2_r1, mcmc.C2_r1, pp.cutoff=0.7)
shiftsumsC2_r2 <- shiftSummaries(chain.C2_r2, mcmc.C2_r2, pp.cutoff=0.7)
shiftsumsC2_r3 <- shiftSummaries(chain.C2_r3, mcmc.C2_r3, pp.cutoff=0.7)
shiftsumsC2_r4 <- shiftSummaries(chain.C2_r4, mcmc.C2_r4, pp.cutoff=0.7)

# Checking the descendents
shiftsumsC2_r1$descendents
shiftsumsC2_r2$descendents
shiftsumsC2_r3$descendents
shiftsumsC2_r4$descendents

# Ploting them
plotShiftSummaries(shiftsumsC2_r1, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsC2_r2, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsC2_r3, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsC2_r4, single.plot=TRUE, label.pts = FALSE)

### Combining the chains that converged. 
### The Gelman plots that show model convergence are located in the code used to run the model in the cluster.
C2_mcmc.chain.list <- list(chain.C2_r1, chain.C2_r2, chain.C2_r4)
C2.mcmc <- combine.chains(chain.list = C2_mcmc.chain.list)

saveRDS(C2.mcmc, "result_mcmc/combined_C2.RDS")

#### TRADITIONAL CLASSIFICATION #####
### In this model, the fixed the locations where shifts could occur and only allowed the intercept to vary. Slope was global

mcmc.T2_r1 = readRDS("result_mcmc/mcmcT2-RUN1.rds")
mcmc.T2_r2 = readRDS("result_mcmc/mcmcT2-RUN2.rds")
mcmc.T2_r3 = readRDS("result_mcmc/mcmcT2-RUN3.rds")
mcmc.T2_r4 = readRDS("result_mcmc/mcmcT2-RUN4.rds")

chain.T2_r1 <- set.burnin(mcmc.T2_r1$load(), 0.3)
chain.T2_r2 <- set.burnin(mcmc.T2_r2$load(), 0.3)
chain.T2_r3 <- set.burnin(mcmc.T2_r3$load(), 0.3)
chain.T2_r4 <- set.burnin(mcmc.T2_r4$load(), 0.3)

### Checking the shifts, just to have an idea if the chains give the same results.
shiftsumsT2_r1 <- shiftSummaries(chain.T2_r1, mcmc.T2_r1, pp.cutoff=0.7)
shiftsumsT2_r2 <- shiftSummaries(chain.T2_r2, mcmc.T2_r2, pp.cutoff=0.7)
shiftsumsT2_r3 <- shiftSummaries(chain.T2_r3, mcmc.T2_r3, pp.cutoff=0.7)
shiftsumsT2_r4 <- shiftSummaries(chain.T2_r4, mcmc.T2_r4, pp.cutoff=0.7)

# Checking the descendents
shiftsumsT2_r1$descendents
shiftsumsT2_r2$descendents
shiftsumsT2_r3$descendents
shiftsumsT2_r4$descendents

# Ploting them
plotShiftSummaries(shiftsumsT2_r1, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsT2_r2, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsT2_r3, single.plot=TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsT2_r4, single.plot=TRUE, label.pts = FALSE)

### Combining the chains that converged. 
### The Gelman plots that show model convergence are located in the code used to run the model in the cluster.
T2_mcmc.chain.list <- list(chain.T2_r1, chain.T2_r3, chain.T2_r4)
T2.mcmc <- combine.chains(chain.list = T2_mcmc.chain.list)

saveRDS(T2.mcmc, "result_mcmc/combined_T2.RDS")

## I'm going to clean the workspace because these things can devour your RAM
rm(list=ls())
 
####### Checking shifts and the differences between the shifts in the models ########
## To do this, we need to load the combined chains, but we also needs the model used to run the mcmc. NOT the individual chain, but the 
## commands and arguments used to build the mcmc chain. That is why I'm also loading one of the mcmcs.

mcmc.11_r1 = readRDS("result_mcmc/mcmc11-RUN1.rds")
chain.11 = readRDS("result_mcmc/combined_11.RDS") 

mcmc.N1_r1 = readRDS("result_mcmc/mcmcN1-RUN1.rds")
chain.N1 = readRDS("result_mcmc/combined_N1.RDS") 

mcmc.NN_r1 = readRDS("result_mcmc/mcmcNN-RUN1.rds")
chain.NN = readRDS("result_mcmc/combined_NN.RDS") 

mcmc.F2_r1 = readRDS("result_mcmc/mcmcF2-RUN1.rds")
chain.F2 = readRDS("result_mcmc/combined_F2.RDS") 

mcmc.C2_r1 = readRDS("result_mcmc/mcmcC2-RUN1.rds")
chain.C2 = readRDS("result_mcmc/combined_C2.RDS") 

mcmc.T2_r1 = readRDS("result_mcmc/mcmcT2-RUN1.rds")
chain.T2 = readRDS("result_mcmc/combined_T2.RDS") 


## Loading the inputed variables for the models, like the tree, the data and body length. We will need that in some cases.
pruned = mcmc.N1_r1$tree
dat = mcmc.N1_r1$dat
eggs = mcmc.N1_r1$pred

### Now, let's build the shift summaries for each model using two posterior probabilit cutoffs 
shiftsumsN1 = shiftSummaries(chain.N1, mcmc.N1_r1, pp.cutoff=0.7)
shiftsumsN1.5 = shiftSummaries(chain.N1, mcmc.N1_r1, pp.cutoff=0.5)

shiftsumsNN = shiftSummaries(chain.NN, mcmc.NN_r1, pp.cutoff=0.7)
shiftsumsNN.5 = shiftSummaries(chain.NN, mcmc.NN_r1, pp.cutoff=0.5)

shiftsumsF2 = shiftSummaries(chain.F2, mcmc.F2_r1, pp.cutoff=0.7)
shiftsumsF2.5 = shiftSummaries(chain.F2, mcmc.F2_r1, pp.cutoff=0.5)

shiftsumsC2 = shiftSummaries(chain.C2, mcmc.C2_r1, pp.cutoff=0.7)
shiftsumsC2.5 = shiftSummaries(chain.C2, mcmc.C2_r1, pp.cutoff=0.5)

shiftsumsT2 = shiftSummaries(chain.T2, mcmc.T2_r1, pp.cutoff=0.7)
shiftsumsT2.5 = shiftSummaries(chain.T2, mcmc.T2_r1, pp.cutoff=0.5)

#ploting shifts
plotShiftSummaries(shiftsumsN1, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsN1.5, lwd=2, single.plot=TRUE, label.pts=FALSE)

plotShiftSummaries(shiftsumsNN, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsNN.5, lwd=2, single.plot=TRUE, label.pts=FALSE)

plotShiftSummaries(shiftsumsF2, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsF2.5, lwd=2, single.plot=TRUE, label.pts=FALSE)

plotShiftSummaries(shiftsumsC2, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsC2.5, lwd=2, single.plot=TRUE, label.pts=FALSE)

plotShiftSummaries(shiftsumsT2, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsT2.5, lwd=2, single.plot=TRUE, label.pts=FALSE)


plotSimmap.mcmc(chain.N1, burnin=0.3, pp.cutoff=0.7)
plotSimmap.mcmc(chain.N1, burnin=0.3, pp.cutoff=0.5)

plotSimmap.mcmc(chain.NN, burnin=0.3, pp.cutoff=0.7)
plotSimmap.mcmc(chain.NN, burnin=0.3, pp.cutoff=0.5)

plotSimmap.mcmc(chain.F2, burnin=0.3, pp.cutoff = 0.7)
plotSimmap.mcmc(chain.F2, burnin=0.3, pp.cutoff = 0.5)

plotSimmap.mcmc(chain.C2, burnin=0.3, pp.cutoff = 0.7)
plotSimmap.mcmc(chain.C2, burnin=0.3, pp.cutoff = 0.5)

plotSimmap.mcmc(chain.T2, burnin=0.3, pp.cutoff = 0.7)
plotSimmap.mcmc(chain.T2, burnin=0.3, pp.cutoff = 0.5)


phenogram.density(pruned,dat,burnin=0.3,chain.11,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.11,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.N1,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.N1,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.NN,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.NN,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.F2,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.F2,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.C2,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.C2,pp.cutoff = 0.5)

phenogram.density(pruned,dat,burnin=0.3,chain.T2,pp.cutoff = 0.7)
phenogram.density(pruned,dat,burnin=0.3,chain.T2,pp.cutoff = 0.5)


#### Now we can check how each parameter shifts across the tree.
## NN is the only model in which the beta varies. Thus, there's no point in plotting the slope for other models.
plotBranchHeatMap(pruned, chain.N1, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))

plotBranchHeatMap(pruned, chain.NN, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))
plotBranchHeatMap(pruned, chain.NN, variable="beta_logBody", burnin=0.3, pal=rainbow, legend_ticks=seq(0.2, 1.5, 0.1))

plotBranchHeatMap(pruned, chain.F2, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))

plotBranchHeatMap(pruned, chain.C2, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))

plotBranchHeatMap(pruned, chain.T2, variable="theta", burnin=0.3, pal=rainbow, legend_ticks=seq(-6,6,1))


#### STEPSTONES ####

### Here, we will load the stepstones andcalculate the bayes factor to know which model is better.

ss.11 = readRDS("result_mcmc/ss_11.rds")
ss.N1 = readRDS("result_mcmc/ss_N1.rds")
ss.NN = readRDS("result_mcmc/ss_NN.rds")
ss.F2 = readRDS("result_mcmc/ss_F2.rds")
ss.T2 = readRDS("result_mcmc/ss_T2.rds")
ss.C2 = readRDS("result_mcmc/ss_C2.rds")


mlnL = c("11"=ss.11$lnr, "N1"=ss.N1$lnr, "NN"=ss.NN$lnr, "M"=ss.F$lnr, 
          "T" = ss.T$lnr,  "C" = ss.C$lnr)
mlnL

summ_BF = as.data.frame(c(2*(mlnL[2]-mlnL[1]),2*(mlnL[2]-mlnL[3]),
                          2*(mlnL[2]-mlnL[4]), 2*(mlnL[2]-mlnL[5]),
                          2*(mlnL[2]-mlnL[6])))
rownames(summ_BF) = c("N1 vs 11", "N1 vs NN","N1 vs F","N1 vs T", "N1 vs C")
colnames(summ_BF) = c("Bayes_Factor")
summ_BF


##### FIGURES FOR THE PAPER ######

#This first one is directly related to the bayes factor.

png("figures/clutchsize-Bayesfactor.png", res = 600, unit = 'mm', h = 160, w = 160)
plot(summ_BF$Bayes_Factor, ylim = c(0,2005), xlab="Models", ylab="Bayes Factor", pch=19,
     xaxt = 'n', cex = 2.2, bty = 'l', las = 1)
axis(side = 1, at = c(1:5), labels = rownames(summ_BF))
abline(h=10, lwd = 2, lty = 2, col = 'red')
dev.off()

## Now, the figures with the shift sums.
# Given N1 was the best model, all figures will be based on it.
mcmc.N1 = readRDS("result_mcmc/combined_N1.RDS")
mcmc.N1_r1 = readRDS("result_mcmc/mcmcN1-RUN1.rds")

pruned = mcmc.N1_r1$tree
dat = mcmc.N1_r1$dat
body = mcmc.N1_r1$pred

plot(dat ~ body)

shiftsums.N1 <- shiftSummaries(mcmc.N1, mcmc.N1_r1, pp.cutoff=0.7)

pdf("figures/simmap.eggnumber-0.7pcutoff.pdf",h=18,w=12)
plotSimmap.mcmc(mcmc.N1, burnin=0.3, edge.type = "regimes", pp.cutoff=0.7, no.margin=TRUE, 
                pal=colorRampPalette(colors=c("#E0B27F","#6B2DE0", "#5AE073", "#0A9424")))
dev.off()

pdf("figures/simmap.eggnumber-0.5pcutoff.pdf",h=18,w=12)
plotSimmap.mcmc(mcmc.N1, burnin=0.3, edge.type = "regimes", pp.cutoff=0.5, no.margin=TRUE, 
                pal=colorRampPalette(colors=c("#E0B27F","#6B2DE0", "#5AE073", "#0A9424")))
dev.off()


shiftsums.N1$descendents
shiftsums.N1$regressions
shiftsums.N1$cladesummaries[[1]]$densities

estimates = shiftsums.N1$regressions

### Loading the data file just to get the correct name of the species.
eggs2 = read.csv("data/LH_data_reorder.csv",h=T,sep=';')
colnames(eggs2)[1] = "species"
eggs2$species = gsub(" ", "_", eggs2$species)

eggs2$logBody = log(eggs2$body.size.avg,10)
eggs2$logFec = log(eggs2$fecundity.avg,10)


#### GRAPH FOR pp 0.7 ####

shiftsums.N1 <- shiftSummaries(mcmc.N1, mcmc.N1_r1, pp.cutoff=0.7)
plotShiftSummaries(shiftsums.N1, single.plot=TRUE, label.pts = FALSE)

shiftsums.N1$descendents

cherax = eggs2[with(eggs2, species %in% shiftsums.N1$descendents[[2]]),]
euastac = eggs2[with(eggs2, species %in% shiftsums.N1$descendents[[3]]),]
astacop = eggs2[with(eggs2, species %in% shiftsums.N1$descendents[[4]]),]

dat = data.frame(eggs2$logFec)

row.names(dat) = eggs2$species
dat = as.matrix(dat)[,1]


logBody = as.matrix(data.frame(eggs2$logBody)) 
rownames(logBody) = eggs2$species

dat = data.frame(eggs2$logFec)

row.names(dat) = eggs2$species
dat = as.matrix(dat)[,1]

### Getting the colors right
#> pal=colorRampPalette(colors=c("#E0B27F","#6B2DE0", "#5AE073", "#0A9424"))
#> pal(7)
#[1] "#E0B27F" "#A56FAF" "#6B2DE0" "#6286A9" "#5AE073" "#32BA4B" "#0A9424"

#png("figures/carapace-clutch-pp0.7.png",res=600,units='mm',h=150,w=200)
tiff("figures/carapace-clutch-pp0.7.tiff",res=600,units='mm',h=150,w=200,
     compression = 'lzw')

plot(dat ~ logBody, cex = 2, las = 1, bty = 'l', 
     ylab = "Clutch size (log10)", xlab = "Body length (log10)",
     pch = 21, bg = "#E0B27F")

points(x = cherax$logBody, y = cherax$logFec, cex = 2, pch = 21, bg = "#0A9424")
points(x = euastac$logBody, y = euastac$logFec, cex = 2, pch = 21, bg =  "#6286A9")
points(x = astacop$logBody, y = astacop$logFec, cex = 2, pch = 21, bg = "#0215a1")

estimates = shiftsums.N1$regressions

curve(estimates[1,1]+(estimates[1,2]*x), add =T, lwd = 3, col = "#807160", from = 1.12, to = 1.85)
curve(estimates[2,1]+(estimates[2,2]*x), add =T, lwd = 3, col = '#0A9424', from = 1.25, to = 1.8)
curve(estimates[3,1]+(estimates[3,2]*x), add =T, lwd = 3, col = "#6286A9", from = 1.65, to = 1.75)
curve(estimates[4,1]+(estimates[4,2]*x), add =T, lwd = 3, col = "#0215a1", from = 2.1, to = 2.2)

dev.off()



#### PP CUTOFF 0.5 #####

shiftsums.N1.05 <- shiftSummaries(mcmc.N1, mcmc.N1_r1, pp.cutoff=0.5)

shiftsums.N1.05$descendents

proc = eggs2[with(eggs2, species %in% shiftsums.N1.05$descendents[[2]]),]
cherax = eggs2[with(eggs2, species %in% shiftsums.N1.05$descendents[[3]]),]
euastac = eggs2[with(eggs2, species %in% shiftsums.N1.05$descendents[[4]]),]
euastac2 = eggs2[with(eggs2, species %in% shiftsums.N1.05$descendents[[5]]),]
euastac3 = eggs2[with(eggs2, species %in% shiftsums.N1.05$descendents[[6]]),]
astacop = eggs2[with(eggs2, species %in% shiftsums.N1.05$descendents[[7]]),]


#"#E0B27F" "#A56FAF" "#6B2DE0" "#6286A9" "#5AE073" "#32BA4B" "#0A9424"

png("figures/carapace-clutch-pp0.5.png",res=600,units='mm',h=150,w=200)

tiff("figures/carapace-clutch-pp0.5.tiff",res=600,units='mm',h=150,w=200,
     compression = 'lzw')
plot(dat ~ logBody, cex = 2, las = 1, bty = 'l', 
     ylab = "Clutch size (log10)", xlab = "Body length (log10)",
     pch = 21, bg = "#E0B27F")

points(x = proc$logBody, y = proc$logFec, cex = 2, pch = 21, bg = "#A56FAF")
points(x = cherax$logBody, y = cherax$logFec, cex = 2, pch = 21, bg = "#0A9424")
points(x = euastac$logBody, y = euastac$logFec, cex = 2, pch = 21, bg =  "#8ef336")
points(x = euastac2$logBody, y = euastac2$logFec, cex = 2, pch = 21, bg =  "#6286A9")
points(x = euastac3$logBody, y = euastac3$logFec, cex = 2, pch = 21, bg =  "#6B2DE0")
points(x = astacop$logBody, y = astacop$logFec, cex = 2, pch = 21, bg = "#0215a1")

estimates = shiftsums.N1.05$regressions

curve(estimates[1,1]+(estimates[1,2]*x), add =T, lwd = 3, col = "#807160", from = 1.12, to = 1.85)
curve(estimates[2,1]+(estimates[2,2]*x), add =T, lwd = 3, col = '#A56FAF', from = 1.25, to = 1.8)
curve(estimates[3,1]+(estimates[3,2]*x), add =T, lwd = 3, col = "#0A9424", from = 1.2, to = 1.8)
curve(estimates[4,1]+(estimates[4,2]*x), add =T, lwd = 3, col = "#8ef336", from = 1.75, to = 1.85)
curve(estimates[5,1]+(estimates[5,2]*x), add =T, lwd = 3, col = "#6286A9", from = 1.68, to = 1.75)
curve(estimates[6,1]+(estimates[6,2]*x), add =T, lwd = 3, col = "#6B2DE0", from = 1.68, to = 1.72)
curve(estimates[7,1]+(estimates[7,2]*x), add =T, lwd = 3, col = "#0215a1", from = 2.1, to = 2.2)
dev.off()