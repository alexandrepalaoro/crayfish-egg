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

ME = sd(eggs$logFec)/sqrt(length(eggs$logFec))

prior.11 <- make.prior(pruned, plot.prior = F, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_logBody="dnorm",
                                  dsb="fixed", dk="fixed", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_logBody=list(mean=1, sd=0.5),
                                  dtheta=list(mean=mean(eggs$logFec), sd=2*sd(eggs$logFec))),
                       fixed = list( k = 0, sb = numeric(0)))

D11 = list(alpha=0.1, sig2=0.1, beta_logBody=1.2, k=0, theta=-0.5, slide=1)

model.11 <- makeBayouModel(dat ~ logBody, rjpars = c(), 
                           tree=pruned, dat=dat, pred=logBody, prior=prior.11, D=D11, SE = ME)

mcmc.11 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.11$model, prior=prior.11, startpar=model.11$startpar, 
                          file.dir="result_mcmc", outname="model11", plot.freq=NULL, ticker.freq = 2000, 
                          samp = 200, SE = ME)

c_chain11 = readRDS("result_mcmc/combined_11.RDS")

Bk <- qbeta(seq(0,1, length.out=2), 0.3, 1)

registerDoParallel(cores = 10)

ss.11 <- mcmc.11$steppingstone(10000, c_chain11, Bk, burnin=0.25, plot=FALSE)

saveRDS(ss.11, "result_mcmc/ss_11.RDS")