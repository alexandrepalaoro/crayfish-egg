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


fixed.loc = list()
fixed.loc$sb = c(193,196,158,115,109 ,101,69,64,51,49)
fixed.loc$loc = c(23.44988, 41.81040,3.58560, 35.33024,12.10403,9.15705,11.65844,13.52917,16.04743,12.36620)
fixed.loc$k = 10
fixed.loc$t2 = c(2,3,4,5,6,7,8,9,10,11)

prior_dists <- list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy"
                    , dbeta_logBody = "dnorm"
                    , dsb = "fixed", dk = "fixed", dtheta = "dnorm")

prior_param <- list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                    dbeta_logBody = list(mean = 1, sd = 0.5),
                    dtheta = list(mean=mean(eggs.na$logsize), sd=2*sd(eggs.na$logsize)))


prior_fixed <- list(k = fixed.loc$k, sb = fixed.loc$sb
                    , loc = fixed.loc$loc, t2 = fixed.loc$t2)

fixPrior = make.prior(tree = pruned, plot.prior = T, dists = prior_dists,
                      param = prior_param, fixed = prior_fixed)


DF = list(alpha = 0.1, sig2 = 0.1, beta_logBody = 1, k = c(1,1), theta = 0, slide = 1)


model.F <- makeBayouModel(dat ~ logBody, rjpars = c("theta")
                          , tree = pruned, dat = dat, pred = logBody
                          , prior = fixPrior, D = DF, SE = ME)

model.F$startpar$loc <- prior_fixed$loc



mcmc.F <- bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.F$model
                         , prior = fixPrior, startpar = model.F$startpar
                         , file.dir="result_mcmc", outname="modelF_r001"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)

mcmc.F_r1 <- bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                            , model = model.F$model
                            , prior = fixPrior, startpar = model.F$startpar
                            , file.dir="result_mcmc", outname="modelF_r1"
                            , plot.freq=NULL, ticker.freq = 2000, samp = 200)
mcmc.F_r2 <- bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                            , model = model.F$model
                            , prior = fixPrior, startpar = model.F$startpar
                            ,file.dir="result_mcmc", outname="modelF_r2"
                            , plot.freq=NULL, ticker.freq = 2000, samp = 200)
mcmc.F_r3 <- bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                            , model = model.F$model
                            , prior = fixPrior, startpar = model.F$startpar
                            ,file.dir="result_mcmc", outname="modelF_r3"
                            , plot.freq=NULL, ticker.freq = 2000, samp = 200)
mcmc.F_r4 <- bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                            , model = model.F$model
                            , prior = fixPrior, startpar = model.F$startpar
                            ,file.dir="result_mcmc", outname="modelF_r4"
                            , plot.freq=NULL, ticker.freq = 2000, samp = 200)


set.seed(sample(1:1000,1))
mcmc.F_r1$run(5000000)
saveRDS(mcmc.F_r1,file="result_mcmc/mcmcF-RUN1.rds")
set.seed(sample(1:1000,1))
mcmc.F_r2$run(5000000)
saveRDS(mcmc.F_r2,file="result_mcmc/mcmcF-RUN2.rds")
set.seed(sample(1:1000,1))
mcmc.F_r3$run(5000000)
saveRDS(mcmc.F_r3,file="result_mcmc/mcmcF-RUN3.rds")
set.seed(sample(1:1000,1))
mcmc.F_r4$run(5000000)
saveRDS(mcmc.F_r4,file="result_mcmc/mcmcF-RUN4.rds")

chain.F_r1 <- set.burnin(mcmc.F_r1$load(), 0.3)
chain.F_r2 <- set.burnin(mcmc.F_r2$load(), 0.3)
chain.F_r3 <- set.burnin(mcmc.F_r3$load(), 0.3)
chain.F_r4 <- set.burnin(mcmc.F_r4$load(), 0.3)

## MODEL F
## PARAMETER: lnL
#PASTIS 
pdf("result_mcmc/Rgelman_lnL_F.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
RlnL <- gelman.R("lnL", chain1=chain.F_r1, chain2=chain.F_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.F_r1, chain2=chain.F_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.F_r1, chain2=chain.F_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.F_r2, chain2=chain.F_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.F_r2, chain2=chain.F_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.F_r3, chain2=chain.F_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()

# PARAMETER: alpha
pdf("result_mcmc/Rgelman_alpha_F.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
Ralpha <- gelman.R("alpha", chain1=chain.F_r1, chain2=chain.F_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.F_r1, chain2=chain.F_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.F_r1, chain2=chain.F_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.F_r2, chain2=chain.F_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.F_r2, chain2=chain.F_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.F_r3, chain2=chain.F_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()

# PARAMETER: sig2
pdf("result_mcmc/Rgelman_sig2_F.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
Rsig2 <- gelman.R("sig2", chain1=chain.F_r1, chain2=chain.F_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.F_r1, chain2=chain.F_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.F_r1, chain2=chain.F_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.F_r2, chain2=chain.F_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.F_r2, chain2=chain.F_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.F_r3, chain2=chain.F_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()