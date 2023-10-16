require(bayou)
library(foreach)
library(doParallel)

eggs = read.csv("data/LH_data_reorder.csv",h=T,sep=';')
colnames(eggs)[1] = "species"
eggs$species = gsub(" ", "_", eggs$species)

eggs$logBody = log(eggs$body.size.avg,10)
eggs$logFec = log(eggs$fecundity.avg,10)

rownames(eggs) = eggs$species

phy = read.tree("data/stern.2017.ml.new.txt")

pruned = drop.tip(phy,phy$tip.label[-match(eggs$species,phy$tip.label)])

name.check(pruned,eggs)

logBody = as.matrix(data.frame(eggs$logBody)) 
rownames(logBody) = eggs$species

dat = data.frame(eggs$logFec)

row.names(dat) = eggs$species
dat = as.matrix(dat)[,1]

fixed.loc = list()
fixed.loc$sb = c(311,242,248,180,133,127,115,85, 76, 70)
fixed.loc$loc = c(63.449220, 6.895874,38.670861, 2.470758,33.685630,11.752591 ,4.745127,18.122865,14.370299, 3.629326)
fixed.loc$k = 10
fixed.loc$t2 = c(2,2,2,2,2,2,2,2,2,2)

ME = sd(eggs$logFec)/sqrt(length(eggs$logFec))

prior_dists <- list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy"
                    , dbeta_logBody = "dnorm"
                    , dsb = "fixed", dk = "fixed", dtheta = "dnorm")

prior_param <- list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                    dbeta_logBody = list(mean = 1, sd = 0.5),
                    dtheta=list(mean=mean(eggs$logFec), sd=2*sd(eggs$logFec)))

prior_fixed <- list(k = fixed.loc$k, sb = fixed.loc$sb
                    , loc = fixed.loc$loc, t2 = fixed.loc$t2)

fixPrior = make.prior(tree = pruned, plot.prior = F, dists = prior_dists,
                      param = prior_param, fixed = prior_fixed)

DF = list(alpha = 0.1, sig2 = 0.1, beta_logBody = 1, k = c(1,1), theta = 0, slide = 1)


model.F2 <- makeBayouModel(dat ~ logBody, rjpars = c("theta")
                           , tree = pruned, dat = dat, pred = logBody
                           , prior = fixPrior, D = DF, SE = ME)
model.F2$startpar$loc <- prior_fixed$loc

mcmc.F2 <- bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.F2$model
                         , prior = fixPrior, startpar = model.F2$startpar
                         , file.dir="result_mcmc", outname="modelF2"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)

c_chainF2 = readRDS("result_mcmc/combined_F2.RDS")

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

registerDoParallel(cores = 10)

ss.F2 <- mcmc.F2$steppingstone(10000000, c_chainF2, Bk, burnin=0.3, plot=FALSE)

saveRDS(ss.F2, "result_mcmc/ss_F2.RDS")