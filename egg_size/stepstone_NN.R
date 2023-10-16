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

prior.NN <- make.prior(pruned, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_logBody="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_logBody=list(mean=1, sd=0.5),
                                  dk=list(lambda=20, kmax=100), 
                                  dtheta=list(mean=mean(eggs.na$logsize), sd=2*sd(eggs.na$logsize)))
)

DNN = list(alpha=1, sig2=1, beta_logBody=1.2, k=c(1,1), theta=0, slide=1)

model.NN <- makeBayouModel(dat ~ logBody, rjpars = c("theta", "logBody"),  
                           tree=pruned, dat=dat, pred=logBody, prior=prior.NN, D=DNN, SE = ME)

mcmc.NN <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.NN$model, prior=prior.NN, startpar=model.NN$startpar, 
                             file.dir="result_mcmc",outname="modelNN", plot.freq=NULL,ticker.freq = 2000, samp = 200, SE = ME)      

c_chainNN = readRDS("result_mcmc/combined_NN.RDS")

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

registerDoParallel(cores = 10)

ss.NN <- mcmc.NN$steppingstone(10000000, c_chainNN, Bk, burnin=0.3, plot=FALSE)

saveRDS(ss.NN, "result_mcmc/ss_NN.RDS")