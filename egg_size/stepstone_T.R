require(bayou)
require(tidyr)
library(foreach)
library(doParallel)

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
fixed.loc$sb = c(197,193,189,182,159,157,149,141,129,119,113,114,109,103,101,92, 69, 61, 51, 49 )
fixed.loc$loc = c(32.030124, 44.627675, 42.768821, 53.506668, 10.056648,  8.952268, 17.657497,  
                  9.939179,  9.415557,  7.614151,  9.559301, 7.927578, 17.678361, 22.349309,
                  15.348641,  5.408437, 20.525066,  2.575059, 27.883294, 26.023330)
fixed.loc$k = 20
fixed.loc$t2 = 2:21

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


DT = list(alpha = 0.1, sig2 = 0.1, beta_logBody = 1, k = c(1,1), theta = 0, slide = 1)


model.T <- makeBayouModel(dat ~ logBody, rjpars = c("theta")
                          , tree = pruned, dat = dat, pred = logBody
                          , prior = fixPrior, D = DT, SE = ME)

model.T$startpar$loc <- prior_fixed$loc



mcmc.T <- bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.T$model
                         , prior = fixPrior, startpar = model.T$startpar
                         , file.dir="result_mcmc", outname="modelT_r001"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)

c_chainT = readRDS("result_mcmc/combined_T.RDS")

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

registerDoParallel(cores = 10)

ss.T <- mcmc.T$steppingstone(10000000, c_chainT, Bk, burnin=0.3, plot=FALSE)

saveRDS(ss.T, "result_mcmc/ss_T.RDS")