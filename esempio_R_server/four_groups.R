library(rjags)
source("src/generate_data.R")

run.jags <- function(obs, nsamples) {
  var = c(
    "fromIdio", "ndpAtomsMean", "ndpComp", "idiosincComp", "sharedAtomsMean",
    "ndpAtomsSigma", "sharedAtomsSigma", "lambdaShared","gamma", "alpha0", "wShared",
    "lambdaNdp", "wNdp", "wNdpTop", "lndpW", "sharedComp", "mean", "std", "meanpred",
    "stdPred", "ypred")
  num_groups = nrow(obs)
  data_ = list(y=obs, H=50, L=30, num_samples=nsamples, numGroups=num_groups)
  mod = jags.model(file="src/lndp.jags", data=data_, n.chains=1, n.adapt=15000)
  
  burnin=5000
  nsamples=5000
  if (burnin > 0)
    update(mod, n.iter=burnin)
  
  coda.samples(
    mod, var, nsamples, thin = 5)
}


ys = matrix(nrow=4, ncol=1000)
ys[1, ] = normal.mixture.samples(c(5, 10), c(0.9, 0.1), 0.6, 100)
ys[2, ] = normal.mixture.samples(c(5, 10), c(0.9, 0.1), 0.6, 100)
ys[3, ] = normal.mixture.samples(c(5, 0), c(0.1, 0.9), 0.6, 100)
ys[4, ] = normal.mixture.samples(c(5, 0), c(0.1, 0.9), 0.6, 100)

num_samples = c(100, 100, 100, 100)

chains = run.jags(ys, num_samples)
saveRDS(chains, file="chains_four_groups.RData")
