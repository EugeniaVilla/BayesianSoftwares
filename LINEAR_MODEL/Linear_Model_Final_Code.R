################################################################

##############    LIBRARIES - PACKAGES  #######################

###############################################################

## LIBRARY
library(coda)
library(tictoc)
library(mvtnorm)

# Jags
library(rjags)

# Stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))

# Nimble 
library(nimble) 

# Plots
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# library(purrr)
# library(ggsci)

# Set Working Directory


###############################################################

######################    DATASET   ###########################

###############################################################

## Set the seed for the Random Number Generator
set.seed(431234)

## Create some simulated data

N=100 #number of observations
K<-4  #length of our parameters 16 or also 120 

# Matrix X of covariates
X <- matrix(nrow=N, ncol=K)
X[,1] <- rep(1,N) #intercept
for(i in 2:K){
  X[,i] <- rnorm(n = N, mean = i^2 ,sd = i*2) #covariates 
}

# for K=120
# X <- matrix(nrow=N, ncol=K)
# X[,1] <- rep(1,N) #intercept
# for(i in 2:K){
#   X[,i] <- rnorm(n = N, mean = 0 ,sd = 1) #covariates 
# }
# X <- abs(X)

# TRUE values
beta <- seq(1,K,1) #[1 2 3 4]
sigmasquare <- 5
tau <- sigmasquare^(-1) #Since Jags uses precision rather than variance

# Y: response
Y <- rep(0,N)
for (i in 1:N){
  Y[i] <- rnorm(n=1, mean=X[i,]%*%beta,sd=sigmasquare^0.5)
}

################ MODEL
# Recall that: LINEAR MODEL
# Y_i=(x_i')*beta+eps_i where eps_i ~ iid N(0, sigmasquared)

# Data
N <- length(Y)
p <- dim(X)[2]

# Initial values
b0 <- rep(2,K)
sigma2 <- 1



################################################################

#################    CONJUGATE - MODEL (ZELLNER)   #############

################################################################
# Recall that: ZELLNER PRIOR for LINEAR MODEL
# beta|sigmasquared ~ N_p(beta0, sigmasquared*B0)
# sigmasquared ~ InvGamma(nu0/2,nu0*sigma0/2)
# where p=dimension of beta, including intercept
#       beta0=rep(0,p)
#       B0=c*(t(X) %*% X)^(-1)

# Costant values for Conjugate Model with Zellner Prior
sigma0 <- 1 
nu0 <- 0.0001 
B0 <- solve(t(X) %*% X) #precomputing B0 to save computational time 
# It's common practice to use 'c' like parameter since its value is unkown in the real situation.
# Then you can rerun each software with the hypothetical value achieved by mcmc list.
Omega <- sigmasquare*B0
# B0=diag(nrow=120) for K=120 case

################### JAGS
# Compilation Model phase
model_C<-
  "model{
  for(i in 1:N){
    for (j in 1:K){
    a[i,j] <- X[i,j]*beta[j] 
    }
    mu[i]<-sum(a[i,]) # deterministic node
    ## likelihood
    Y[i] ~ dnorm(mu[i],tau)  # stochastic node
  }

  ## priors
  #REMEMBER Jags uses precision matrix
  sigmasq <- pow(tau,-1) # deterministic node
  OmegaInverse <- inverse(c*Omega)
  beta[1:K] ~ dmnorm(rep(0,K),OmegaInverse[1:K,1:K]) # stochastic node
  tau ~ dgamma(tau0/2,tau0*sigma0/2) # stochastic node
  c ~ dunif(0,1000000) # stochastic node
}"
model0<-textConnection(model_C)

## Compilation phase
{tic()
  JAGS0<-jags.model(file=model0,
                    data=list('N'=N,'K'= p, 'X'=X,'Y'=Y,'tau0'= nu0, 'sigma0' = sigma0, 'Omega' = Omega),
                    inits=list('beta'= beta,'tau'= 1/sigma2), #Let software choose the inital value for parameter 'c'
                    n.chains=4,
                    n.adapt=10000) #n.adapth=burn-in
  ## Sampling phase
  ## Run with CODA
  CODA0<-coda.samples(model=JAGS0,
                      variable.names=c('beta','tau', 'sigmasq','c'), 
                      n.iter=50000,
                      thin=10)
  toc()}

## Funtions for manipulating jags.model.objects
list.samplers(JAGS0)
coef(JAGS0)

## Results by CODA
#head(CODA0)
summary(CODA0)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA0)

## AutoCorrelation Funtion
x11()
acfplot(CODA0)


## Effective Sample Size and Rhat
effectiveSize(CODA0)#ness
gelman.diag(CODA0, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)#Rhat

dev.off()
dev.off()

################### STAN

inits <- function() 
{
  list(
    beta = b0, 
    sigma2=sigma2)
} 

data_lm_zellner <-list(N = N, 
                       y=Y,
                       p = p, 
                       X=as.matrix(X),
                       B0=B0,
                       b0=rep(0,p),
                       sigma0=sigma0,
                       nu0=nu0
) 

{tic()
  LM_Z1 <- stan(file = "LinearModel_Zellner_cparam.stan", 
                data = data_lm_zellner,
                chains = 4, 
                iter = 100000, 
                warmup = 50000, 
                thin= 10, 
                seed = 42, 
                init = inits,
                algorithm = 'NUTS')
  toc()}

save(LM_Z1, file="LM_Z1.dat")

print(LM_Z1, pars = c('beta','sigma2','c'))

mcmcModelZ1<-As.mcmc.list(LM_Z1)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModelZ1[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModelZ1[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModelZ1[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModelZ1[[1]][,4]) # traceplot of beta[4]
x11()
coda::traceplot(mcmcModelZ1[[1]][,5]) # traceplot of sigma2

## AutoCorrelation Funtion
x11()
acf(mcmcModelZ1[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModelZ1[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModelZ1[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModelZ1[[1]][,4]) #autocorrelation for beta[4]
x11()
acf(mcmcModelZ1[[1]][,5]) #autocorrelation for sigma2

## Effective Sample Size and Rhat
effectiveSize(mcmcModelZ1)
gelman.diag(mcmcModelZ1, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

graphics.off() #shuts down all open graphics devices



################### NIMBLE

linearCodeZ <- nimbleCode({
  for(i in 1:n){
    for(j in 1:k){
      a[i,j] <- X[i,j]*beta[j]
    }
    mu[i] <- sum(a[i,])
    # likelihood
    Y[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## priors
  for(i in 1:k){    # with NIMBLE, if we want to define a vector inside the model, we have to do a for cicle!
    b0[i] <- 0
  }
  #B0[1:k,1:k] <- c*(inverse(t(X[1:n,1:k])%*%X[1:n,1:k])) 
  tau0 <- 0.001
  sigma0 <- 1
  c ~ dunif(0,1000000)
  sigmasq ~ dinvgamma(tau0*0.5,tau0*sigma0*0.5)
  tau <- pow(sigmasq,-1)  ## residual std dev
  sigma <- pow(sigmasq,0.5)
  cova[1:k,1:k] <- sigmasq*c*B0[1:k,1:k]
  beta[1:k] ~ dmnorm(b0[1:k],cov=cova[1:k,1:k])   # we cannot pass an expression to dmnorm so we have to create z before!
  
})

linearConsts <- list(n =N,k=p)
linearInits <-  list(beta=b0,sigmasq=sigma2,tau = 1/sigma2)
linearData <- list(X=X,Y=Y,B0=B0)


{tic()
  linear <- nimbleModel(code = linearCodeZ,dimensions = list(a=c(N,p)),
                        name = "linear",constants = linearConsts, data = linearData, inits = linearInits)
  Clinear <- compileNimble(linear)
  linearConf <- configureMCMC(linear, print = TRUE)
  linearConf$addMonitors(c("tau", "sigmasq", "beta","c"))
  linearMCMC <- buildMCMC(linearConf)
  ClinearMCMC <- compileNimble(linearMCMC, project = linear)
  {tic()
    Psamples <- runMCMC(ClinearMCMC, niter=100000, nburnin=50000, thin=10,
                        nchains=4,samplesAsCodaMCMC = TRUE,summary=TRUE)
    toc()}
  toc()}


Psamples$summary

## Diagnostic by CODA
## Traceplot and Density plot
x11()
traceplot(Psamples[["samples"]][["chain1"]][,1], xlab = "iteration", ylab = expression(alpha),main="Beta1 traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,5],xlab = "iteration", ylab = expression(sbeta),main="sigmasq traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,6],xlab = "iteration", ylab = expression(sbeta),main="tau traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,7],xlab = "iteration", ylab = expression(sbeta),main="c traceplot")

## AutoCorrelation Funtion
x11() 
acf(Psamples[["samples"]][["chain1"]][,1],main="Beta1 autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,5],main="sigmasq autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,6],main="tau autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,7],main="c autocorrelation")


## Effective Sample Size and Rhat
effectiveSize(Psamples[[1]])
gelman.diag(Psamples[[1]], confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

graphics.off()


###############################################################

#################    CONJUGATE HIERARCHICAL - MODEL   #########

###############################################################
# beta|sigmasquared ~ N_p(beta0, Sigma)
# Sigma ~ InvWishart(S0^(-1),eta0)
# sigmasquared ~ InvGamma(nu0/2,nu0*sigma0/2) 
# where p=p, dimension of beta, including intercept
#       beta0=rep(0,p)
#       S0=c*((t(X) %*% X)^(-1))
#       eta0=p+2, degree of freedom (eta0>p+1)

# Costant values
nu0 <- 0.0001
sigma0 <- 1
c <- 100
B0 <- c * solve(t(X) %*% X) #B0=diag(120)*c
Omega <- sigmasquare*B0 #S0^(-1)
Sigma_beta <- diag(rep(1,p))

########################## Jags
# Compilation Model phase
model_C<-
  "model{
  for(i in 1:N){
    for (j in 1:K){
    a[i,j] <- X[i,j]*beta[j] 
    }
    mu[i]<-sum(a[i,]) # deterministic node
    ## likelihood
    Y[i] ~ dnorm(mu[i],tau)  # stochastic node
  }

  ## priors
  #REMEMBER Jags uses precision matrix
  sigmasq <- pow(tau,-1) # deterministic node
  beta[1:K] ~ dmnorm(rep(0,K),Sigma[1:K,1:K]) # stochastic node
  Sigma ~ dwish(OmegaInverse,K+2) # stochastic node
  tau ~ dgamma(tau0/2,tau0*sigma0/2) # stochastic node
}"
model1<-textConnection(model_C)

## Compilation phase
{tic()
  JAGS1<-jags.model(file=model1,
                    data=list('N'=N,'K'= p, 'X'=X,'Y'=Y,'tau0'= nu0, 'sigma0' = sigma0, 'OmegaInverse' = B0/sigmasquare),
                    inits=list('beta'= b0,'tau'= 1/sigma2,'Sigma'=Sigma_beta),
                    n.chains=4,
                    n.adapt=25000) #n.adapth=burn-in
  ## Sampling phase
  ## Run with CODA
  
  CODA1<-coda.samples(model=JAGS1,
                      variable.names=c('beta','tau', 'sigmasq'), 
                      n.iter=25000,
                      thin=10)
  toc()}

## Funtions for manipulating jags.model.objects
list.samplers(JAGS1)
coef(JAGS1)

## Results by CODA
#head(CODA1)
summary(CODA1)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA1)

## AutoCorrelation Funtion
x11()
acfplot(CODA1)


## Effective Sample Size and Rhat
effectiveSize(CODA1)#ness
gelman.diag(CODA1, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)#Rhat

dev.off()
dev.off()

############################ Stan

data_lm_wishart <-list(N = N, 
                       y = Y,
                       p = p, 
                       X=as.matrix(X),
                       B0=B0,
                       b0=rep(0,p),
                       sigma0=sigma0,
                       nu0=nu0,
                       Omega =Omega
) 

inits_wishart <- function() 
{
  list(
    beta = b0, 
    sigma2=sigma2,
    Sigma_beta =Sigma_beta
  )
}

## run stan model
{tic()
  LM <- stan(file = "LinearModel_Wishart.stan", 
             data = data_lm_wishart,
             chains = 4, 
             iter = 50000, 
             warmup = 25000,  
             thin= 10, 
             seed = 42, 
             init = inits_wishart,
             control = list(adapt_delta = 0.95),
             algorithm = 'NUTS')
  toc()}

save(LM, file="LM.dat")

print(LM, pars = c('beta','sigma2'))

mcmcModel<-As.mcmc.list(LM)
is.mcmc.list(mcmcModel)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel[[1]][,4]) # traceplot of beta[4]
x11()
coda::traceplot(mcmcModel[[1]][,5]) # traceplot of sigma2


## AutoCorrelation Funtion
x11()
acf(mcmcModel[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel[[1]][,4]) #autocorrelation for beta[4]
x11()
acf(mcmcModel[[1]][,5]) #autocorrelation for sigma2

## Effective Sample Size and Rhat
effectiveSize(mcmcModel[[1]][,c(1:5,122)])# nummber of effective sample size
gelman.diag(mcmcModel, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)#Rhat
#potential scale reduction factor

graphics.off()

######################### Nimble 
linearCode <- nimbleCode({
  for(i in 1:N){
    for(j in 1:K){
      a[i,j] <- X[i,j]*beta[j]
    }
    mu[i] <- sum(a[i,])
    # likelihood
    Y[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## priors
  beta[1:K] ~ dmnorm(beta0[1:K],cov=Sigma[1:K,1:K])
  Sigma[1:K,1:K] ~ dinvwish(S = Omega[1:K,1:K], df = K+2) #Omega is defined outside
  for(i in 1:K){    #with NIMBLE, if we want to define a vector inside the model, we have to do a for cicle!
    beta0[i] <- 2 # beta0 = [0 0 0 0]
  }
  sigma <- pow(sigmasq,0.5)
  sigmasq ~ dinvgamma(nu0,nu0*sigma0/2)
})

linearConsts <- list(c=c,N=N,K=K,nu0=nu0,sigma0=sigma0)
linearInits <- list( beta=b0, sigmasq=sigma2, Sigma=Sigma_beta)
linearData <- list(X=X, Y=Y, Omega=Omega)

{tic()
  linear <- nimbleModel(code = linearCode,dimensions = list(a=c(N,p)),  # I have to specify the dimension !
                        name = "linear", constants = linearConsts, data = linearData, inits = linearInits)
  Clinear <- compileNimble(linear)
  linearConf <- configureMCMC(linear, print = TRUE)
  linearConf$addMonitors(c( "sigmasq", "beta","Sigma"))
  
  linearMCMC <- buildMCMC(linearConf)
  ClinearMCMC <- compileNimble(linearMCMC, project = linear)
  {tic()
    Psamples <- runMCMC(ClinearMCMC, niter=50000, nburnin=25000, thin=10, nchains=4,samplesAsCodaMCMC = TRUE,summary=TRUE)
    toc()}
  toc()}


Psamples$summary

## Effective Sample Size and Rhat
effectiveSize(Psamples[[1]])# nummber of effective sample size
gelman.diag(Psamples[[1]], confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)#Rhat

## Diagnostic by CODA
## Traceplot and Density plot
x11()
traceplot(Psamples[["samples"]][["chain1"]][,1], xlab = "iteration", ylab = expression(alpha),main="Beta1 traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,5],xlab = "iteration", ylab = expression(sbeta),main="sigmasq traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,6],xlab = "iteration", ylab = expression(sbeta),main="tau traceplot")

## AutoCorrelation Funtion
x11()   
acf(Psamples[["samples"]][["chain1"]][,1],main="Beta1 autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,5],main="sigmasq autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,6],main="tau autocorrelation")

graphics.off()



##########################################################

##############   NOT CONJUGATE - MODEL  #########################

##########################################################
# Recall that: NOT CONJUGATE PRIOR for LINEAR MODEL
# beta ~ N_p(0,100) # improper prior
# sigmasquared ~ Unif(0,10000) # improper prior
# where p=p, dimension of beta, including intercept

########################## Jags
# Compilation Not Conjugate Model phase
model_NC<-
  "model{
  for(i in 1:N){
    for (j in 1:K){
    a[i,j] <- X[i,j]*beta[j] 
    }
    mu[i]<-sum(a[i,]) # deterministic node
    ## likelihood
    Y[i] ~ dnorm(mu[i],tau)  # stochastic node
  }
  
  ## priors
  for(j in 1:K){
    beta[j] ~ dnorm(0,0.01) # stochastic node
    # remember that JAGS wants the precision
  }  
  sigmasq ~ dunif(0,10000) # stochastic node
  tau <- pow(sigmasq,-1) # deterministic node   
}"
model2<-textConnection(model_NC)

## Compilation phase
{tic()
  JAGS2<-jags.model(file=model2,
                    data=list('N'=N,'K'= p, 'X'=X,'Y'=Y),
                    inits=list('beta'=b0,'sigmasq'= sigma2),
                    n.chains=4,
                    n.adapt=50000) #n.adapth=burn-in
  
  ## Sampling phase
  ## Run with CODA
  CODA2<-coda.samples(model=JAGS2,
                      variable.names=c('beta','tau', 'sigmasq'), 
                      n.iter=50000,
                      thin=10)
  toc()}

## Funtions for manipulating jags.model.objects
list.samplers(JAGS2)
coef(JAGS2)

## Results by CODA
#head(CODA2)
summary(CODA2)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA2)

## AutoCorrelation Funtion
x11()
acfplot(CODA2)


## Effective Sample Size and Rhat
effectiveSize(CODA2)#ness
gelman.diag(CODA2, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)#Rhat


dev.off()
dev.off()

######################### Stan
data_lm_nc <-list(N = N, 
                  y = Y,
                  p = p, 
                  X = as.matrix(X)
) 

# Initial values
inits_nc <- function() 
{
  list(
    beta = b0, 
    sigma=sigma2)
}

{tic()
  LM_NC <- stan(file = "LinearModel_NotConjugate.stan", 
                data = data_lm_nc,
                chains = 4, 
                iter = 100000, 
                warmup = 50000, 
                thin= 10, 
                seed = 42, 
                init = inits_nc,
                algorithm = 'NUTS')
  toc()}

save(LM_NC, file="LM_NC.dat")

print(LM_NC, pars = c('beta','sigma2'))

mcmcModel_NC<-As.mcmc.list(LM_NC)

## Results by CODA
summary(mcmcModel_NC)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel_NC[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_NC[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_NC[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_NC[[1]][,4]) # traceplot of beta[4]
x11()
coda::traceplot(mcmcModel_NC[[1]][,5]) # traceplot of sigma2

## AutoCorrelation Funtion
x11()
acf(mcmcModel_NC[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_NC[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_NC[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_NC[[1]][,4]) #autocorrelation for beta[4]
x11()
acf(mcmcModel_NC[[1]][,5]) #autocorrelation for sigma2

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_NC)#ness
gelman.diag(mcmcModel_NC, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)#Rhat

graphics.off()

######################### Nimble 
linearCode <- nimbleCode({
  for(i in 1:n){
    for(j in 1:k){
      a[i,j] <- X[i,j]*beta[j]
    }
    mu[i] <- sum(a[i,])
    # likelihood
    Y[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## priors
  for(j in 1:k){
    #beta[j] ~ dnorm(0,sd=10)
    beta[j] ~ dflat()
  }
  sigmasq ~ dflat()   # dflat() <- improper prior
  #sigmasq ~ dunif(0,10000)
  sigma <- pow(sigmasq,0.5)
})

linearConsts <- list(n=N,k=p)
linearInits <-  list(beta=b0,sigmasq=sigma2)
linearData <- list(X=X,Y=Y)

{tic()
  linear <- nimbleModel(code = linearCode,dimensions = list(a=c(100,4)),  # I have to specify the dimension !
                        name = "linear",constants = linearConsts, data = linearData, inits = linearInits)
  Clinear <- compileNimble(linear)
  linearConf <- configureMCMC(linear, print = TRUE)
  linearConf$addMonitors(c("sigmasq", "beta"))
  
  linearMCMC <- buildMCMC(linearConf)
  ClinearMCMC <- compileNimble(linearMCMC, project = linear)
  {tic()
    Psamples <- runMCMC(ClinearMCMC, niter=100000, nburnin=50000, thin=10, nchains=4,samplesAsCodaMCMC = TRUE,summary=TRUE)
    toc()}
  toc()}


Psamples$summary

## Effective Sample Size and Rhat
gelman.diag(Psamples[[1]], confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
effectiveSize(Psamples[[1]])     # effective sample size

## Diagnostic by CODA
## Traceplot and Density plot
x11()
traceplot(Psamples[["samples"]][["chain1"]][,1], xlab = "iteration", ylab = expression(alpha),main="Beta1 traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,5],xlab = "iteration", ylab = expression(sbeta),main="sigmasq traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,6],xlab = "iteration", ylab = expression(sbeta),main="tau traceplot")

## AutoCorrelation Funtion
x11()
acf(Psamples[["samples"]][["chain1"]][,1],main="Beta1 autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,5],main="sigmasq autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,6],main="tau autocorrelation")

graphics.off()
