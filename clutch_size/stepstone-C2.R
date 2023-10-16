require(bayou)
library(foreach)
library(doParallel)

eggs = read.csv("data/LH_data_reorder.csv",h=T,sep=';')
colnames(eggs)[1] = "species"
eggs$species = gsub(" ", "_", eggs$species)

eggs$logBody = scale(log(eggs$body.size.avg,10),scale=F)
eggs$logFec = scale(log(eggs$fecundity.avg,10),scale=F)

rownames(eggs) = eggs$species

phy = read.tree("data/stern.2017.ml.new.txt")

pruned = drop.tip(phy,phy$tip.label[-match(eggs$species,phy$tip.label)])

name.check(pruned,eggs)

logBody = as.matrix(data.frame(eggs$logBody)) 
rownames(logBody) = eggs$species

dat = data.frame(eggs$logFec)

row.names(dat) = eggs$species
dat = as.matrix(dat)[,1]

ME = sd(eggs$logFec)/sqrt(length(eggs$logFec))

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

prior_dists <- list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy"
                    , dbeta_logBody = "dnorm"
                    , dsb = "fixed", dk = "fixed", dtheta = "dnorm")

prior_param <- list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                    dbeta_logBody = list(mean = 1, sd = 0.5),
                    dtheta = list(mean=mean(eggs$logFec), sd=2*sd(eggs$logFec)))

prior_fixed <- list(k = fixed.loc$k, sb = fixed.loc$sb
                    , loc = fixed.loc$loc, t2 = fixed.loc$t2)

fixPrior = make.prior(tree = pruned, plot.prior = T, dists = prior_dists,
                      param = prior_param, fixed = prior_fixed)

DC2 = list(alpha = 0.1, sig2 = 0.1, beta_logBody = 1, k = c(1,1), theta = 0, slide = 1)


model.C2 <- makeBayouModel(dat ~ logBody, rjpars = c("theta")
                          , tree = pruned, dat = dat, pred = logBody
                          , prior = fixPrior, D = DC2, SE = ME)

model.C2$startpar$loc <- prior_fixed$loc


mcmc.C2 <- bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.C2$model
                         , prior = fixPrior, startpar = model.C2$startpar
                         , file.dir="result_mcmc", outname="modelC2_r001"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)

c_chainC2 = readRDS("result_mcmc/combined_C2.RDS")

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)


registerDoParallel(cores = 10)

ss.C2 <- mcmc.C2$steppingstone(10000000, c_chainC2, Bk, burnin=0.25, plot=FALSE)

saveRDS(ss.C2, "result_mcmc/ss_C2.RDS")