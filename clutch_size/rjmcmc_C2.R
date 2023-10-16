require(bayou)

eggs = read.csv("data/LH_data_reorder.csv",h=T,sep=';')
colnames(eggs)[1] = "species"
eggs$species = gsub(" ", "_", eggs$species)

eggs$logBody = log(eggs$body.size.avg)
eggs$logFec = log(eggs$fecundity.avg)

rownames(eggs) = eggs$species

phy = read.tree("data/stern.2017.ml.new.txt")

pruned = drop.tip(phy,phy$tip.label[-match(eggs$species,phy$tip.label)])

name.check(pruned,eggs)

logBody = as.matrix(data.frame(eggs$logBody)) 
rownames(logBody) = eggs$species

dat = data.frame(eggs$logFec)

row.names(dat) = eggs$species
dat = as.matrix(dat)[,1]

### Fixing the location of shifts

fixed.loc = list()
fixed.loc$sb = c(311, 305, 304, 287, 285, 281, 279, 273, 271, 267, 257, 254,
                 245, 234, 238, 248, 183, 181, 161, 147, 137, 133, 127, 121,
                 119, 110,  85,  73,  70, 58)
fixed.loc$loc = c(61.7176408,  7.1203635,5.9869824, 62.4581189, 50.4352449, 52.9961288,
                  45.1854001, 30.0412998,7.1592884, 47.0239855,26.9870879,  4.2750957,
                  65.0832588, 26.5003134,27.6491938, 34.3505650,25.2332220, 11.5563953,
                  9.2256726, 14.2655274,6.6489389, 27.2547237,22.3981506, 26.6351231,
                  19.5476628,  3.3756591,22.1723523,  2.6166516,4.3067221, 13.6881458)
fixed.loc$k = 30
fixed.loc$t2 = c(rep(2,29),1)


ME = sd(eggs$logFec)/sqrt(length(eggs$logFec))

prior_dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy"
                    , dbeta_logBody = "dnorm"
                    , dsb = "fixed", dk = "fixed", dtheta = "dnorm")

prior_param = list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                    dbeta_logBody = list(mean = 1, sd = 0.5),
                    dtheta = list(mean=mean(eggs$logFec), sd=2*sd(eggs$logFec)))

prior_fixed = list(k = fixed.loc$k, sb = fixed.loc$sb
                    , loc = fixed.loc$loc, t2 = fixed.loc$t2)

fixPrior = make.prior(tree = pruned, plot.prior = T, dists = prior_dists,
                      param = prior_param, fixed = prior_fixed)

DC2 = list(alpha = 0.1, sig2 = 0.1, beta_logBody = 1, k = c(1,1), theta = 0, slide = 1)


model.C2 = makeBayouModel(dat ~ logBody, rjpars = c("theta")
                          , tree = pruned, dat = dat, pred = logBody
                          , prior = fixPrior, D = DC2, SE = ME)

model.C2$startpar$loc = prior_fixed$loc

mcmc.C2 = bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.C2$model
                         , prior = fixPrior, startpar = model.C2$startpar
                         , file.dir="result_mcmc", outname="modelC2_r001"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)

mcmc.C2_r1 = bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.C2$model
                         , prior = fixPrior, startpar = model.C2$startpar
                         , file.dir="result_mcmc", outname="modelC2_r1"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)
mcmc.C2_r2 = bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.C2$model
                         , prior = fixPrior, startpar = model.C2$startpar
                         ,file.dir="result_mcmc", outname="modelC2_r2"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)
mcmc.C2_r3 = bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.C2$model
                         , prior = fixPrior, startpar = model.C2$startpar
                         ,file.dir="result_mcmc", outname="modelC2_r3"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)
mcmc.C2_r4 = bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.C2$model
                         , prior = fixPrior, startpar = model.C2$startpar
                         ,file.dir="result_mcmc", outname="modelC2_r4"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)


set.seed(sample(1:1000,1))
mcmc.C2_r1$run(5000000)
saveRDS(mcmc.C2_r1,file="result_mcmc/mcmcC2-RUN1.rds")
set.seed(sample(1:1000,1))
mcmc.C2_r2$run(5000000)
saveRDS(mcmc.C2_r2,file="result_mcmc/mcmcC2-RUN2.rds")
set.seed(sample(1:1000,1))
mcmc.C2_r3$run(5000000)
saveRDS(mcmc.C2_r3,file="result_mcmc/mcmcC2-RUN3.rds")
set.seed(sample(1:1000,1))
mcmc.C2_r4$run(5000000)
saveRDS(mcmc.C2_r4,file="result_mcmc/mcmcC2-RUN4.rds")


chain.C2_r1 = set.burnin(mcmc.C2_r1$load(), 0.3)
chain.C2_r2 = set.burnin(mcmc.C2_r2$load(), 0.3)
chain.C2_r3 = set.burnin(mcmc.C2_r3$load(), 0.3)
chain.C2_r4 = set.burnin(mcmc.C2_r4$load(), 0.3)

## MODEL F
## PARAMETER: lnL
#PASTIS 
pdf("result_mcmc/Rgelman_lnL_C2.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
RlnL <- gelman.R("lnL", chain1=chain.C2_r1, chain2=chain.C2_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.C2_r1, chain2=chain.C2_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.C2_r1, chain2=chain.C2_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.C2_r2, chain2=chain.C2_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.C2_r2, chain2=chain.C2_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
RlnL <- gelman.R("lnL", chain1=chain.C2_r3, chain2=chain.C2_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()

# PARAMETER: alpha
pdf("result_mcmc/Rgelman_alpha_C2.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
Ralpha <- gelman.R("alpha", chain1=chain.C2_r1, chain2=chain.C2_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.C2_r1, chain2=chain.C2_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.C2_r1, chain2=chain.C2_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.C2_r2, chain2=chain.C2_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.C2_r2, chain2=chain.C2_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain.C2_r3, chain2=chain.C2_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()

# PARAMETER: sig2
pdf("result_mcmc/Rgelman_sig2_C2.pdf", height = 11, width = 8.5)
par(mfrow=c(2,3))
Rsig2 <- gelman.R("sig2", chain1=chain.C2_r1, chain2=chain.C2_r2, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.C2_r1, chain2=chain.C2_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.C2_r1, chain2=chain.C2_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.C2_r2, chain2=chain.C2_r3, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.C2_r2, chain2=chain.C2_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain.C2_r3, chain2=chain.C2_r4, plot=TRUE, type="n", ylim=c(0.9, 2))
dev.off()

