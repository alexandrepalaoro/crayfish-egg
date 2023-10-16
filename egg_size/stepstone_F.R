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

c_chainF = readRDS("result_mcmc/combined_F.RDS")

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

registerDoParallel(cores = 10)

ss.F <- mcmc.F$steppingstone(10000000, c_chainF, Bk, burnin=0.3, plot=FALSE)

saveRDS(ss.F, "result_mcmc/ss_F.RDS")