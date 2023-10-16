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

prior.N1 <- make.prior(pruned, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_logBody="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_logBody=list(mean=1, sd=0.5),
                                  dk=list(lambda=20, kmax=100),
                                  dtheta=list(mean=mean(eggs.na$logsize), sd=2*sd(eggs.na$logsize)))
)


DN1 = list(alpha=0.1, sig2=0.1, beta_logBody=1, k=1, theta=0, slide=1)

model.N1 <- makeBayouModel(dat ~ logBody, rjpars = c("theta"),  
                           tree=pruned, dat=dat, pred=logBody, prior=prior.N1, D=DN1, SE = ME)

mcmc.N1_r1 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, 
                          file.dir="result_mcmc",outname="modelN1_r1", plot.freq=NULL,ticker.freq = 2000, samp = 200, SE = ME)      
mcmc.N1_r2 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, 
                          file.dir="result_mcmc",outname="modelN1_r2", plot.freq=NULL,ticker.freq = 2000, samp = 200, SE = ME)
mcmc.N1_r3 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, 
                          file.dir="result_mcmc",outname="modelN1_r3", plot.freq=NULL,ticker.freq = 2000, samp = 200, SE = ME)
mcmc.N1_r4 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, 
                          file.dir="result_mcmc",outname="modelN1_r4", plot.freq=NULL,ticker.freq = 2000, samp = 200, SE = ME)

set.seed(sample(1:1000,1))
mcmc.N1_r1$run(5000000)
saveRDS(mcmc.N1_r1,file="result_mcmc/mcmcN1-RUN1.rds")
set.seed(sample(1:1000,1))
mcmc.N1_r2$run(5000000)
saveRDS(mcmc.N1_r2,file="result_mcmc/mcmcN1-RUN2.rds")
set.seed(sample(1:1000,1))
mcmc.N1_r3$run(5000000)
saveRDS(mcmc.N1_r3,file="result_mcmc/mcmcN1-RUN3.rds")
set.seed(sample(1:1000,1))
mcmc.N1_r4$run(5000000)
saveRDS(mcmc.N1_r4,file="result_mcmc/mcmcN1-RUN4.rds")

chain.N1_r1 <- set.burnin(mcmc.N1_r1$load(), 0.3)
chain.N1_r2 <- set.burnin(mcmc.N1_r2$load(), 0.3)
chain.N1_r3 <- set.burnin(mcmc.N1_r3$load(), 0.3)
chain.N1_r4 <- set.burnin(mcmc.N1_r4$load(), 0.3)

## MODEL N1
## PARAMETER: lnL
#PASTIS 
pdf("result_mcmc/Rgelman_lnL_N1.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
RlnL <- gelman.R("lnL", chain1=chain.N1_r1, chain2=chain.N1_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.N1_r1, chain2=chain.N1_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.N1_r1, chain2=chain.N1_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.N1_r2, chain2=chain.N1_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.N1_r2, chain2=chain.N1_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.N1_r3, chain2=chain.N1_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()

# PARAMETER: alpha
pdf("result_mcmc/Rgelman_alpha_N1.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
Ralpha <- gelman.R("alpha", chain1=chain.N1_r1, chain2=chain.N1_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.N1_r1, chain2=chain.N1_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.N1_r1, chain2=chain.N1_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.N1_r2, chain2=chain.N1_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.N1_r2, chain2=chain.N1_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.N1_r3, chain2=chain.N1_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()

# PARAMETER: sig2
pdf("result_mcmc/Rgelman_sig2_N1.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
Rsig2 <- gelman.R("sig2", chain1=chain.N1_r1, chain2=chain.N1_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.N1_r1, chain2=chain.N1_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.N1_r1, chain2=chain.N1_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.N1_r2, chain2=chain.N1_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.N1_r2, chain2=chain.N1_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.N1_r3, chain2=chain.N1_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()