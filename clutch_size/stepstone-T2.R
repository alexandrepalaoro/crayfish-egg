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
fixed.loc$sb = c(311, 305, 304, 287, 285, 281, 279, 273, 271, 267, 257, 245, 239, 234, 238, 248, 223, 220, 210, 183, 
                 181, 180, 171, 161, 147, 137, 131, 132, 127, 110,  85,  70 , 58)
fixed.loc$loc = c(65.705091,11.556162, 11.217769, 50.135377, 41.098502, 50.273806, 40.943154, 26.902663,  6.227868, 
                  45.487785, 27.286663, 64.801230, 33.498537, 14.159623, 29.958546, 32.252453, 47.504553,  7.752675, 
                  54.454662, 27.985812, 11.730326,  2.078817, 13.282972, 14.221378, 14.166737,  8.459463 , 8.396882,  
                  5.184776,  9.144059,  5.855528, 25.439978,  4.753375, 14.296341)
fixed.loc$k = 33
fixed.loc$t2 = c(2,2,3,3,3,3,3,3,2,2,3,2,3,2,2,2,3,3,3,3,3,2,3,3,3,3,2,3,3,2,3,2,1)

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

DT2 = list(alpha = 0.1, sig2 = 0.1, beta_logBody = 1, k = c(1,1), theta = 0, slide = 1)


model.T2 <- makeBayouModel(dat ~ logBody, rjpars = c("theta")
                          , tree = pruned, dat = dat, pred = logBody
                          , prior = fixPrior, D = DT2, SE = ME)

model.T2$startpar$loc <- prior_fixed$loc


mcmc.T2 <- bayou.makeMCMC(tree = pruned, dat = dat, SE = ME, pred = logBody
                         , model = model.T2$model
                         , prior = fixPrior, startpar = model.T2$startpar
                         , file.dir="result_mcmc", outname="modelT2_r001"
                         , plot.freq=NULL, ticker.freq = 2000, samp = 200)

c_chainT2 = readRDS("result_mcmc/combined_T2.RDS")

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

registerDoParallel(cores = 10)

ss.T2 <- mcmc.T2$steppingstone(10000000, c_chainT2, Bk, burnin=0.25, plot=FALSE)

saveRDS(ss.T2, "result_mcmc/ss_T2.RDS")