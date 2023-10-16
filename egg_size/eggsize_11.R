require(bayou)
require(tidyr)

eggs = read.csv("data/LH_data_reorder.csv",h=T,sep=';')
colnames(eggs)[1] = "species"
eggs$species = gsub(" ", "_", eggs$species)

eggs$logBody = log(eggs$body.size.avg,10)
eggs$logsize = log(eggs$egg.diam.avg,10)

eggs.na = eggs %>% drop_na(logsize)

rownames(eggs.na) = eggs.na$species

phy = read.tree("data/stern.2017.ml.new.txt")

pruned = drop.tip(phy,phy$tip.label[-match(eggs.na$species,phy$tip.label)])

name.check(pruned,eggs.na)

logBody = as.matrix(data.frame(eggs.na$logBody)) 
rownames(logBody) = eggs.na$species

dat = data.frame(eggs.na$logsize)

row.names(dat) = eggs.na$species
dat = as.matrix(dat)[,1]

ME = sd(eggs.na$logsize)/sqrt(length(eggs.na$logsize))


prior.11 <- make.prior(pruned, plot.prior = F, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_logBody="dnorm",
                                  dsb="fixed", dk="fixed", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_logBody=list(mean=1, sd=0.5),
                                  dtheta=list(mean=mean(eggs.na$logsize), sd=2*sd(eggs.na$logsize))),
                       fixed = list( k = 0, sb = numeric(0)))

D11 = list(alpha=0.1, sig2=0.1, beta_logBody=1, k=0, theta=0.3, slide=1)

model.11 <- makeBayouModel(dat ~ logBody, rjpars = c(), 
                           tree=pruned, dat=dat, pred=logBody, prior=prior.11, D=D11, SE = ME)

mcmc.11_r1 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.11$model, prior=prior.11, startpar=model.11$startpar, 
                          file.dir="result_mcmc", outname="model11_r1", plot.freq=NULL, ticker.freq = 2000, 
                          samp = 200, SE = ME)
mcmc.11_r2 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.11$model, prior=prior.11, startpar=model.11$startpar, 
                          file.dir="result_mcmc", outname="model11_r2", plot.freq=NULL, ticker.freq = 2000, 
                          samp = 200, SE = ME)
mcmc.11_r3 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.11$model, prior=prior.11, startpar=model.11$startpar, 
                          file.dir="result_mcmc", outname="model11_r3", plot.freq=NULL, ticker.freq = 2000, 
                          samp = 200, SE = ME)
mcmc.11_r4 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.11$model, prior=prior.11, startpar=model.11$startpar, 
                          file.dir="result_mcmc", outname="model11_r4", plot.freq=NULL, ticker.freq = 2000, 
                          samp = 200, SE = ME)

set.seed(sample(1:1000,1))
mcmc.11_r1$run(5000000)
saveRDS(mcmc.11_r1,file="result_mcmc/mcmc11-RUN1.rds")
set.seed(sample(1:1000,1))
mcmc.11_r2$run(5000000)
saveRDS(mcmc.11_r2,file="result_mcmc/mcmc11-RUN2.rds")
set.seed(sample(1:1000,1))
mcmc.11_r3$run(5000000)
saveRDS(mcmc.11_r3,file="result_mcmc/mcmc11-RUN3.rds")
set.seed(sample(1:1000,1))
mcmc.11_r4$run(5000000)
saveRDS(mcmc.11_r4,file="result_mcmc/mcmc11-RUN4.rds")

chain.11_r1 <- set.burnin(mcmc.11_r1$load(), 0.3)
chain.11_r2 <- set.burnin(mcmc.11_r2$load(), 0.3)
chain.11_r3 <- set.burnin(mcmc.11_r3$load(), 0.3)
chain.11_r4 <- set.burnin(mcmc.11_r4$load(), 0.3)


## MODEL 11
## PARAMETER: lnL
#PASTIS 
pdf("result_mcmc/Rgelman_lnL_11.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
RlnL <- gelman.R("lnL", chain1=chain.11_r1, chain2=chain.11_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.11_r1, chain2=chain.11_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.11_r1, chain2=chain.11_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.11_r2, chain2=chain.11_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.11_r2, chain2=chain.11_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.11_r3, chain2=chain.11_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()

# PARAMETER: alpha
pdf("result_mcmc/Rgelman_alpha_11.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
Ralpha <- gelman.R("alpha", chain1=chain.11_r1, chain2=chain.11_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.11_r1, chain2=chain.11_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.11_r1, chain2=chain.11_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.11_r2, chain2=chain.11_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.11_r2, chain2=chain.11_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.11_r3, chain2=chain.11_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()

# PARAMETER: sig2
pdf("result_mcmc/Rgelman_sig2_11.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
Rsig2 <- gelman.R("sig2", chain1=chain.11_r1, chain2=chain.11_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.11_r1, chain2=chain.11_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.11_r1, chain2=chain.11_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.11_r2, chain2=chain.11_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.11_r2, chain2=chain.11_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.11_r3, chain2=chain.11_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()