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

prior.N1 <- make.prior(pruned, plot.prior = F, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_logBody="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_logBody=list(mean=1, sd=0.5),
                                  dk=list(lambda=20, kmax=100),
                                  dtheta=list(mean=mean(eggs$logFec), sd=2*sd(eggs$logFec)))
)


DN1 = list(alpha=0.1, sig2=0.1, beta_logBody=1, k=1, theta=0, slide=1)

model.N1 <- makeBayouModel(dat ~ logBody, rjpars = c("theta"),  
                           tree=pruned, dat=dat, pred=logBody, prior=prior.N1, D=DN1, SE = ME)

mcmc.N1 <- bayou.makeMCMC(pruned, dat, pred=logBody, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, 
                         file.dir="result_mcmc", outname="modelN1_r001", plot.freq=NULL, ticker.freq = 2000, samp = 200, SE = ME)

c_chainN1 = readRDS("result_mcmc/combined_N1.RDS")

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

registerDoParallel(cores = 10)

ss.N1 <- mcmc.N1$steppingstone(10000000, c_chainN1, Bk, burnin=0.25, plot=FALSE)

saveRDS(ss.N1, "result_mcmc/ss_N1.RDS")