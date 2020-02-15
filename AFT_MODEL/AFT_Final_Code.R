#########################################################################

####################### LIBRARY & PACKAGES ##############################

#########################################################################

library(stats)

#Packages for JAGS and CODA
library(rjags)
library(coda)
library(tictoc)
library(nimble)

#Packages for STAN
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))

# Set Working Directory


#########################################################################

############################## DATASET ##################################

#########################################################################

######### DATA GENERATION FOR AFT MODEL

########################

set.seed(423156)

N <- 1000 # number of elements in the dataset

# simulate covariates 

p<-4
X<-matrix(nrow=N, ncol=p)
X[,1]<-rep(1,N) # intercept
for(i in 2:p){
  X[,i] <- rnorm(n = N, mean = p^0.1 , sd =p/8 ) # covariates, which have negative values
}

# true parameters

beta <- seq(0.2,0.2+(p-1)*0.75, 0.75) # [0.20 0.95 1.70 2.45]
mu <- X %*% beta

sigma <- 1

alpha <- 1/sigma # parameter shape for weibull distribution
lambda <- exp(mu)/((log(2))^(1/alpha)) # parameter scale for weibull distribution

# RECALL:
# The Weibull distribution with shape parameter a 
# and scale parameter b has density given by
# 
# f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a) for x > 0.
# 
# The cumulative distribution function is F(x) = 1 - exp(- (x/b)^a) on x > 0, 
# the mean is E(X) = b G(1 + 1/a), and the Var(X) = b^2 * (G(1 + 2/a) - (G(1 + 1/a))^2).

# simulate censoring and survival times
survt = rweibull(N, shape=alpha, scale = lambda) 
cent = rweibull(N, shape=alpha, scale = lambda)

## observed data:
#censoring indicator
delta <- survt < cent # delta is equal to 0 if data is censured
survt[delta==0] <- cent[delta==0] # censor survival time.

# count number of missing/censored survival times
n_miss <- N-sum(delta)

# data for censored subjects
y_m=survt[delta==0]
X_m=X[delta==0,]

# data for uncensored subjects
y_o=survt[delta==1]
X_o=X[delta==1,]
N_m = n_miss
N_o = N - n_miss

#####################################################################

#################### NOT HIERARCHICAL PRIOR #########################

#####################################################################

# initial values of parameters
beta0=rep(0.1, p) #[0.1 0.1 0.1 0.1]
alpha0=0.1

################# JAGS ########################
# Compilation AFT Model
model_AFT1<-"model {
 for(i in 1:N){
   for(s in 1:p){
   a[i,s] <- X[i,s]*(beta[s])
   }
   mu[i] <- sum(a[i,])
   lambda[i] <- log(2)*exp(-mu[i]*alpha)
   # likelihood
   censured[i] ~ dinterval(t[i], cent[i])
   t[i] ~ dweib(alpha, lambda[i])
   }
   #Priors
   for(i in 1:p){
   beta[i] ~ dnorm(0,0.1)
   }
   alpha ~ dexp(1)
}"
model1<-textConnection(model_AFT1)

## Compilation phase
{tic()
  JAGS1<-jags.model(file=model1,
                    data=list('N'=N,'X'=X,'p'= p,'cent'=cent, 't'=survt),
                    inits=list('beta'=beta0,'alpha'=alpha0),
                    n.chains=2,
                    n.adapt=5000) #n.adapth=burn-in
  ## Sampling phase
  ## Run with CODA
  CODA1<-coda.samples(model=JAGS1,
                      variable.names=c('beta','alpha'), 
                      n.iter=5000,
                      thin=5) 
  toc()}

## Funtions for manipulating jags.model.objects
list.samplers(JAGS1)
coef(JAGS1)

## Results by CODA
summary(CODA1)

alpha 
beta
## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA1)

## AutoCorrelation Funtion
x11()
acfplot(CODA1)


## Effective Sample Size
effectiveSize(CODA1)

dev.off()
dev.off()
################# STAN ########################

data_aft <-list(p=p,
                N_m = N_m,
                X_m=as.matrix(X_m),
                y_m=y_m,
                N_o = N_o,
                X_o=as.matrix(X_o),
                y_o=y_o)


#initialization of the parameters

inits <- function() 
{
  list(
    beta = beta0, 
    alpha=alpha0
  )
}

# run stan model
{tic()
AFT <- stan(file = "AFT.stan", 
            data = data_aft,
            chains = 2, 
            iter = 50000, 
            warmup = 10000, 
            thin= 5, 
            seed = 42, 
            init = inits,
            algorithm = 'NUTS')
toc()}

save(AFT, file="AFT.dat")

print(AFT, pars=c('beta','alpha'))

mcmcModel_AFT<-As.mcmc.list(AFT)


## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel_AFT[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_AFT[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_AFT[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_AFT[[1]][,4]) # traceplot of beta[4]
x11()
coda::traceplot(mcmcModel_AFT[[1]][,5]) # traceplot of alfa


## AutoCorrelation Funtion
x11()
acf(mcmcModel_AFT[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_AFT[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_AFT[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_AFT[[1]][,4]) #autocorrelation for beta[4]
x11()
acf(mcmcModel_AFT[[1]][,5]) #autocorrelation for alfa

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_AFT)
gelman.diag(mcmcModel_AFT, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

##################### NIMBLE #######################

################# NOT-HIERARCHICAL


# initial values of parameters
beta0=rep(0.1, p) #[0.1 0.1 0.1 0.1]
alpha0=0.1

AFT_NH_Code <- nimbleCode({
   #Likelihood     
   for(i in 1:N){
      for(s in 1:p){
         a[i,s] <- X[i,s]*(beta[s])
      }
      mu[i] <- sum(a[i,])
      lambda[i] <- log(2)*exp(-mu[i]*alpha)
      censured[i] ~ dinterval(t[i], cent[i])
      t[i] ~ dweib(alpha, lambda[i])
   }
   #Priors
   for(i in 1:p){
      beta[i] ~ dnorm(0,sd=sqrt(10))
   }
   alpha ~ dexp(1)
})


aftConsts <- list(N=N,p=p)
aftInits <-  list(beta=beta0,alpha=alpha0)
aftData <- list(X=X,cent=cent,t=survt)

{tic()
   aftNH <- nimbleModel(code = AFT_NH_Code,dimensions = list(a = c(N,p)),
                        name = "aftNH",constants = aftConsts, data = aftData, inits = aftInits)
   Caft <- compileNimble(aftNH)
   aftConf <- configureMCMC(aftNH, print = TRUE)
   aftConf$addMonitors(c("beta","alpha"))
   aftMCMC <- buildMCMC(aftConf)
   CaftMCMC <- compileNimble(aftMCMC, project = aftNH)
   {tic()
      Psamples <- runMCMC(CaftMCMC, niter=10000, nburnin=5000, thin=10,
                          nchains=2,samplesAsCodaMCMC = TRUE,summary=TRUE)
      toc()}
   toc()}

Psamples$summary

gelman.diag(Psamples[[1]], confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
effectiveSize(Psamples[[1]])

########## Traceplot
x11()
traceplot(Psamples[["samples"]][["chain1"]][,1], xlab = "iteration", ylab = expression(alpha),main="Beta1 traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,5],xlab = "iteration", ylab = expression(sbeta),main="sigmasq traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,6],xlab = "iteration", ylab = expression(sbeta),main="tau traceplot")

x11()   #-----------   AUTOCORRELATION
acf(Psamples[["samples"]][["chain1"]][,1],main="Beta1 autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,5],main="sigmasq autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,6],main="tau autocorrelation")



graphics.off() #shuts down all open graphics devices

#####################################################################

######################## HIERARCHICAL PRIOR #########################

#####################################################################

#### hyperparameters
mu_hat=1
nu_hat=4
sigma2_hat=5
c=2
d=2

# initial values of parameters
beta0=rep(0.1, p)
alpha0=0.1
mu0=0.5
sigma20=1


################# JAGS ########################
# Compilation AFT Model
model_AFT2<-"model {
 for(i in 1:N){
   for(s in 1:p){
   a[i,s] <- X[i,s]*(beta[s])
   }
   mu[i] <- sum(a[i,])
   lambda[i] <- log(2)*exp(-mu[i]*alpha) 
   # likelihood
   censured[i] ~ dinterval(t[i], cent[i])
   t[i] ~ dweib(alpha, lambda[i])
   }
   #Priors
   for(i in 1:p){
   beta[i] ~ dnorm(mub,tau)
   }
   mub ~ dnorm(mu_0,tau)
   sigma2 <- pow(tau,-1)
   tau ~ dgamma(nu_0/2,nu_0*sigma2_0/2)
   alpha ~ dgamma(c,d)
   
}"
model2<-textConnection(model_AFT2)

## Compilation phase
{tic()
  JAGS2<-jags.model(file=model2,
                    data=list('N'=N,'X'=X,'p'= p,'cent'=cent, 't'=survt,'mu_0'=mu_hat,'sigma2_0'=sigma2_hat,'nu_0'=nu_hat,
                              'c'=c,'d'=d),
                    inits=list('beta'=beta0,'alpha'=alpha0,'tau'=1/sigma20,'mub'=mu0),
                    n.chains=2,
                    n.adapt=5000) #n.adapth=burn-in
  ## Sampling phase
  ## Run with CODA
  CODA2<-coda.samples(model=JAGS2,
                      variable.names=c('beta','alpha','mub','sigma2'), 
                      n.iter=20000,
                      thin=5) 
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


## Effective Sample Size
effectiveSize(CODA2)

dev.off()
dev.off()
################# STAN ########################
data_aft_H <-list(p=p,
                N_m = N_m,
                X_m=as.matrix(X_m),
                y_m=y_m,
                N_o = N_o,
                X_o=as.matrix(X_o),
                y_o=y_o,
                mu_0=mu_hat,
                nu_0=nu_hat,
                sigma2_0=sigma2_hat,
                c=c,
                d=d)


#initialization of the parameters

inits_H <- function() 
{
  list(
    beta = beta0, 
    alpha=alpha0,
    mu=mu0,
    sigma2=sigma20
  )
}

# run stan model
{tic()
  AFT_H <- stan(file = "AFT_Hierarchical.stan", 
              data = data_aft_H,
              chains = 2, 
              iter = 50000, 
              warmup = 10000, 
              thin= 5, 
              seed = 42, 
              init = inits_H,
              algorithm = 'NUTS')
  toc()}

save(AFT_H, file="AFT_H.dat")

print(AFT_H, pars=c('beta','alpha'))

mcmcModel_AFT_H<-As.mcmc.list(AFT_H)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel_AFT_H[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_AFT_H[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_AFT_H[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_AFT_H[[1]][,4]) # traceplot of beta[4]
x11()
coda::traceplot(mcmcModel_AFT_H[[1]][,5]) # traceplot of alfa
x11()
coda::traceplot(mcmcModel_AFT_H[[1]][,6]) # traceplot of mu
x11()
coda::traceplot(mcmcModel_AFT_H[[1]][,7]) # traceplot of sigma2


## AutoCorrelation Funtion
x11()
acf(mcmcModel_AFT_H[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_AFT_H[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_AFT_H[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_AFT_H[[1]][,4]) #autocorrelation for beta[4]
x11()
acf(mcmcModel_AFT_H[[1]][,5]) #autocorrelation for alfa
x11()
acf(mcmcModel_AFT_H[[1]][,6]) #autocorrelation for mu
x11()
acf(mcmcModel_AFT_H[[1]][,7]) #autocorrelation for sigma2

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_AFT_H)
gelman.diag(mcmcModel_AFT_H, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

##################### NIMBLE #######################

################# HIERARCHICAL

#### hyperparameters
mu_hat=1
nu_hat=4
sigma2_hat=5
c=2
d=2

# initial values of parameters
beta0=rep(0.1, p)
alpha0=0.1
mu0=0.5
sigma20=1

AFT_H_Code <- nimbleCode({
   #Likelihood     
   for(i in 1:N){
      for(s in 1:p){
         a[i,s] <- X[i,s]*(beta[s])
      }
      mu[i] <- sum(a[i,])
      lambda[i] <- log(2)*exp(-mu[i]*alpha)
      censured[i] ~ dinterval(t[i], cent[i])
      t[i] ~ dweib(alpha, lambda[i])
   }
   #Priors
   for(i in 1:p){
      beta[i] ~ dnorm(mub,sd=sigma)
   }
   mub ~ dnorm(mu_0,sd=sigma)
   sigma2 ~ dinvgamma(nu_0/2,nu_0*sigma2_0/2)
   sigma <- pow(sigma2,0.5)
   alpha ~ dinvgamma(c,d)
})


aftConsts <- list(N=N,p=p)
aftInits <-  list(beta=beta0,alpha=alpha0,mub=mu0)
aftData <- list(X=X,cent=cent,t=survt,mu_0=mu_hat,sigma2_0=sigma2_hat,nu_0=nu_hat,c=c,d=d)

{tic()
   aftH <- nimbleModel(code = AFT_H_Code,dimensions = list(a = c(N,p)),
                       name = "aftH",constants = aftConsts, data = aftData, inits = aftInits)
   Caft <- compileNimble(aftH)
   aftConf <- configureMCMC(aftH, print = TRUE)
   aftConf$addMonitors(c("beta","alpha","mub","sigma2"))
   aftMCMC <- buildMCMC(aftConf)
   CaftMCMC <- compileNimble(aftMCMC, project = aftH)
   {tic()
      Psamples <- runMCMC(CaftMCMC, niter=10000, nburnin=5000, thin=10,
                          nchains=2,samplesAsCodaMCMC = TRUE,summary=TRUE)
      toc()}
   toc()}

Psamples$summary

gelman.diag(Psamples[[1]], confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
effectiveSize(Psamples[[1]])

########## Traceplot
x11()
traceplot(Psamples[["samples"]][["chain1"]][,1], xlab = "iteration", ylab = expression(alpha),main="Beta1 traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,5],xlab = "iteration", ylab = expression(sbeta),main="sigmasq traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,6],xlab = "iteration", ylab = expression(sbeta),main="tau traceplot")

x11()   #-----------   AUTOCORRELATION
acf(Psamples[["samples"]][["chain1"]][,1],main="Beta1 autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,5],main="sigmasq autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,6],main="tau autocorrelation")


graphics.off() #shuts down all open graphics devices

#####################################################################

#################### NOT INFORMATIVE PRIOR ##########################

#####################################################################

################# JAGS ########################
# Compilation AFT Model
model_AFT2<-"model {
 for(i in 1:N){
   for(s in 1:p){
   a[i,s] <- X[i,s]*(beta[s])
   }
   mu[i] <- sum(a[i,])
   lambda[i] <- log(2)*exp(-mu[i]*alpha) 
   # likelihood
   censured[i] ~ dinterval(t[i], cent[i])
   t[i] ~ dweib(alpha, lambda[i])
   }
   #Priors
   for(i in 1:p){
   beta[i] ~ dunif(-100,100)
   }
   alpha ~ dunif(0.01,100)
}"
model3<-textConnection(model_AFT3)

## Compilation phase
{tic()
  JAGS3<-jags.model(file=model3,
                    data=list('N'=N,'X'=X,'p'= p,'cent'=cent, 't'=survt),
                    inits=list('beta'=beta0,'alpha'=alpha0),
                    n.chains=2,
                    n.adapt=5000) #n.adapth=burn-in
  ## Sampling phase
  ## Run with CODA
  CODA3<-coda.samples(model=JAGS3,
                      variable.names=c('beta','alpha'), 
                      n.iter=20000,
                      thin=5) 
  toc()}

## Funtions for manipulating jags.model.objects
list.samplers(JAGS3)
coef(JAGS3)

## Results by CODA
#head(CODA3)
summary(CODA3)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA3)

## AutoCorrelation Funtion
x11()
acfplot(CODA3)


## Effective Sample Size
effectiveSize(CODA3)

dev.off()
dev.off()
################# STAN ########################

data_aft_NI <-list(p=p,
                N_m = N_m,
                X_m=as.matrix(X_m),
                y_m=y_m,
                N_o = N_o,
                X_o=as.matrix(X_o),
                y_o=y_o)


#initialization of the parameters

inits_NI <- function() 
{
  list(
    beta = beta0, 
    alpha=alpha0
  )
}

# run stan model
{tic()
  AFT_NI <- stan(file = "AFT_Non_Informative.stan", 
              data = data_aft_NI,
              chains = 2, 
              iter = 50000, 
              warmup = 10000, 
              thin= 5, 
              seed = 42, 
              init = inits_NI,
              algorithm = 'NUTS')
  toc()}

save(AFT_NI, file="AFT_NI.dat")

print(AFT_NI, pars=c('beta','alpha'))

mcmcModel_AFT_NI<-As.mcmc.list(AFT_NI)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel_AFT_NI[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_AFT_NI[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_AFT_NI[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_AFT_NI[[1]][,4]) # traceplot of beta[4]
x11()
coda::traceplot(mcmcModel_AFT_NI[[1]][,5]) # traceplot of alfa


## AutoCorrelation Funtion
x11()
acf(mcmcModel_AFT_NI[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_AFT_NI[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_AFT_NI[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_AFT_NI[[1]][,4]) #autocorrelation for beta[4]
x11()
acf(mcmcModel_AFT_NI[[1]][,5]) #autocorrelation for alfa

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_AFT_NI)
gelman.diag(mcmcModel_AFT_NI, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

##################### NIMBLE #######################

################# NOT-INFORMATIVE


AFT_NI_Code <- nimbleCode({
   for(i in 1:N){
      for(s in 1:p){
         a[i,s] <- X[i,s]*(beta[s])
      }
      mu[i] <- sum(a[i,])
      lambda[i] <- log(2)*exp(-mu[i]*alpha) 
      # likelihood
      censured[i] ~ dinterval(t[i], cent[i])
      t[i] ~ dweib(alpha, lambda[i])
   }
   #Priors
   for(i in 1:p){
      beta[i] ~ dunif(-100,100)
   }
   alpha ~ dunif(0.01,100)
})


aftConsts <- list(N=N,p=p)
aftInits <-  list(beta=beta0,alpha=alpha0)
aftData <- list(X=X,cent=cent,t=survt)

{tic()
   aftNI <- nimbleModel(code = AFT_NI_Code,dimensions = list(a = c(N,p)),
                        name = "aftNI",constants = aftConsts, data = aftData, inits = aftInits)
   Caft <- compileNimble(aftNI)
   aftConf <- configureMCMC(aftNI, print = TRUE)
   aftConf$addMonitors(c("beta","alpha"))
   aftMCMC <- buildMCMC(aftConf)
   CaftMCMC <- compileNimble(aftMCMC, project = aftNI)
   {tic()
      Psamples <- runMCMC(CaftMCMC, niter=10000, nburnin=5000, thin=10,
                          nchains=2,samplesAsCodaMCMC = TRUE,summary=TRUE)
      toc()}
   toc()}

Psamples$summary

gelman.diag(Psamples[[1]], confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
effectiveSize(Psamples[[1]])

########## Traceplot
x11()
traceplot(Psamples[["samples"]][["chain1"]][,1], xlab = "iteration", ylab = expression(alpha),main="Beta1 traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,5],xlab = "iteration", ylab = expression(sbeta),main="sigmasq traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,6],xlab = "iteration", ylab = expression(sbeta),main="tau traceplot")

x11()   #-----------   AUTOCORRELATION
acf(Psamples[["samples"]][["chain1"]][,1],main="Beta1 autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,5],main="sigmasq autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,6],main="tau autocorrelation")





graphics.off() #shuts down all open graphics devices
