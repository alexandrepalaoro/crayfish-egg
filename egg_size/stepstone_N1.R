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

mcmc.N1 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, 
                             file.dir="result_mcmc",outname="modelN1", plot.freq=NULL,ticker.freq = 2000, samp = 200, SE = ME)

c_chainN1 = readRDS("result_mcmc/combined_N1.RDS")

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

registerDoParallel(cores = 10)

ss.N1 <- mcmc.N1$steppingstone(10000000, c_chainN1, Bk, burnin=0.25, plot=FALSE)

saveRDS(ss.N1, "result_mcmc/ss_N1.RDS")