#########################################################################

####################### LIBRARY & PACKAGES ##############################

#########################################################################

library(stats)

#Packages for JAGS and CODA
library(rjags)
library(coda)
library(tictoc)

#Packages for STAN
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))

# Set Working Directory

#########################################################################

############################## DATASET ##################################

#########################################################################

######### DATA GENERATION POISSON MODEL
set.seed(215679) 

N <- 100 # number of observations
#N <-1000
p_fix <- 4 # number of prameters
#p_fix <- 16

X <- matrix(nrow=N, ncol=p_fix) # Matrix X of covariates
X[,1] <- rep(1,N) # intercept
for(i in 2:p_fix){
  X[,i] <- rnorm(n = N, mean = p_fix/10 ,sd =p_fix/100 ) #covariates
}
# only case with p_fix=16
# X <- matrix(nrow=N, ncol=p_fix) # Matrix X of covariates
# X[,1] <- rep(1,N) # intercept
# for(i in 2:p_fix){
#   X[,i] <- rnorm(n = N, mean = 0 ,sd =0.1 ) #covariates
# }
# X <- abs(X)  only case with p_fix=16


## parameters
beta <- seq(1,p_fix,1) #[1 2 3 4]
sigmasquare <- rep(4, p_fix)
tau <- rep(1/4, p_fix)

lambda <- rep(0,N)

for (i in 1:N){
  lambda[i] <- exp(X[i,]%*%t(t(beta)))
}
#lambda

## Data
Y <- rep(NA,N)
for (s in 1:N)
{
  Y[s] <- rpois(1, lambda = lambda[s]);  
} 
#Y

## Constant values
beta0<-rep(0,p_fix) 


#########################################################################

####################### NOT HIERARCHICAL PRIOR ##########################

#########################################################################
## Initial values
b_0 <- rep(2,p_fix)

# Compilation Not Hierarchical Model
model_Not_Hierarchical<-
  "model{
  for(i in 1:N){  
    for(j in 1:K){
    a[i,j] <- X[i,j]*beta[j] # deterministic node
    }
    zeta[i] <- sum(a[i,])
    ## likelihood
    Y[i] ~ dpois(exp(zeta[i]))  # stochastic node
  }
  
  ## priors
  for(i in 1:K){
  beta[i] ~ dnorm(beta0[i],tau[i]) # stochastic node
  # remember that JAGS wants the inverse of the precision matrix
  sigmasquared[i]<- pow(tau[i],-1) # deterministic node
  } 
}"
model1<-textConnection(model_Not_Hierarchical)

## Compilation phase
{tic()
  JAGS_NH<-jags.model(file=model1,
                   data=list('N'=N, 'X'=X,'K'= p_fix,'Y'=Y,'beta0'=beta0,'tau'= tau),
                   inits=list('beta'=b_0),
                   n.chains=4,
                   n.adapt=5000) #n.adapth=burn-in #(N=10000 -> idem)
  
  ## Sampling phase
  ## Run with CODA
  CODA_NH<-coda.samples(model=JAGS_NH,
                     variable.names=c('beta'), 
                     n.iter=20000,                  #(N=10000 -> idem)
                     thin=10)
  toc()}

## Funtions for manipulating jags.model.objects
list.samplers(JAGS_NH)
coef(JAGS_NH)

## Results by CODA
#head(CODA)
summary(CODA_NH)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA_NH)

## AutoCorrelation Funtion
x11()
acfplot(CODA_NH)


## Effective Sabmple Size
effectiveSize(CODA_NH)
gelman.diag(CODA_NH)

dev.off()
dev.off()

################# STAN

data_glm <-list(N = N, 
                p_fix=p_fix,
                Y=as.vector(Y),
                X=as.matrix(X),
                sigma2_beta=sigmasquare,
                beta0=beta0
) 


#initialization of the parameters

inits <- function() 
{
  list(
    beta = b_0
  )
}

# run stan model
{tic()
GLM <- stan(file = "GLM.stan", 
            data = data_glm,
            chains = 4, 
            iter = 25000, 
            warmup = 5000, 
            thin= 10, 
            seed = 42, 
            init = inits,
            algorithm = 'NUTS')
toc()}

save(GLM, file="GLM.dat")

print(GLM, pars=('beta'))


mcmcModel_GLM<-As.mcmc.list(GLM)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel_GLM[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_GLM[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_GLM[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_GLM[[1]][,4]) # traceplot of beta[4]


## AutoCorrelation Funtion
x11()
acf(mcmcModel_GLM[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_GLM[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_GLM[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_GLM[[1]][,4]) #autocorrelation for beta[4]

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_GLM)[1:4]
gelman.diag(mcmcModel_GLM, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

graphics.off() #shuts down all open graphics devices

##################### NIMBLE #######################à

model_Not_Hierarchical <- nimbleCode({
  for(i in 1:N){
    for(j in 1:K){
      a[i,j] <- X[i,j]*beta[j]
    }
    zeta[i] <- sum(a[i,]) # deterministic nodes
    ## likelihood
    Y[i] ~ dpois(exp(zeta[i]))  # stochastic nodes
  }
  
  ## priors
  for(i in 1:K){
    beta[i] ~ dnorm(beta0[i],sd=sigma[i])
    sigma[i] <- pow(sigmasquared[i],0.5)
  } 
})


GLMConsts <- list(N=N,K=p_fix)
GLMData <-  list(beta0=beta0,X=X,Y=Y,sigmasquared=sigmasquared)

GLMInits <- list(beta=b_0)

{tic()
  GLM <- nimbleModel(code = model_Not_Hierarchical,dimensions = list(a=c(N,p_fix)),
                     name = "GLM",constants = GLMConsts, data = GLMData, inits = GLMInits)
  C_GLM <- compileNimble(GLM)
  GLMConf <- configureMCMC(GLM, print = TRUE)
  GLMConf$addMonitors(c("beta"))
  
  GLMMCMC <- buildMCMC(GLMConf)
  C_GLMMCMC <- compileNimble(GLMMCMC, project = GLM)
  {tic()
    Psamples <- runMCMC(C_GLMMCMC, niter=10000, nburnin=5000, thin=10, nchains=4,samplesAsCodaMCMC = TRUE,summary=TRUE)
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


#########################################################################

########################### HIERARCHICAL PRIOR ##########################

#########################################################################

## Constant values
a0 <- 1/3
b0 <- 1/8

## Initial values
sigma2_0 <- rep(10, p_fix) #real value [4 4 4 4]
tau_0 <- rep(0.1, p_fix)
b_0 <- rep(2,p_fix) #real value [1 2 3 4]
alpha_0 <- c(5,5) #real value [3 8]


# Compilation Hierarchical Model
model_Hierarchical<-
  "model{
  for(i in 1:N){  
    for(j in 1:K){
    a[i,j] <- X[i,j]*beta[j] # deterministic node
    }
    zeta[i] <- sum(a[i,])
    ## likelihood
    Y[i] ~ dpois(exp(zeta[i]))  # stochastic node
  }
  
  ## priors
  for(i in 1:K){
  beta[i] ~ dnorm(beta0[i],tau[i]) # stochastic node
  # remember that JAGS wants the inverse of the precision matrix
  sigmasquared[i]<- pow(tau[i],-1) # deterministic node
  tau[i] ~ dgamma(alpha[1],alpha[2]) # stochastic node
  }
  alpha[1] ~ dexp(a0) # stochastic node
  alpha[2] ~ dexp(b0) # stochastic node
}"
model2<-textConnection(model_Hierarchical)

## Compilation phase
{tic()

JAGS_H<-jags.model(file=model2,
                 data=list('N'=N, 'X'=X,'K'= p_fix,'Y'=Y,'a0'=a0, 'b0'=b0,'beta0'=beta0),
                 inits=list('beta'=b_0,'tau'=tau_0,'alpha'=alpha_0),
                 n.chains=4,
                 n.adapt=20000) #n.adapth=burn-in

## Sampling phase
## Run with CODA
CODA_H<-coda.samples(model=JAGS_H,
                   variable.names=c('beta', 'sigmasquared','alpha'), 
                   n.iter=20000,
                   thin=10) 
toc()}

## Funtions for manipulating jags.model.objects
list.samplers(JAGS_H)
coef(JAGS_H)

## Results by CODA
#head(CODA_H)
summary(CODA_H)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA_H)

## AutoCorrelation Funtion
x11()
acfplot(CODA_H)


## Effective Sample Size
effectiveSize(CODA_H)
gelman.diag(CODA_H)

dev.off()
dev.off()

##################################
############### STAN


data_glm_h <-list(N = N, 
                p_fix=p_fix,
                Y=as.vector(Y),
                X=as.matrix(X),
                beta0=beta0,
                a0=a0,
                b0=b0
) 


#initialization of the parameters


inits_h <- function() 
{
  list(
    beta = b_0, 
    sigma2_beta=sigma2_0,
    a=alpha_0[1],
    b=alpha_0[2]
  )
}

# run stan model
{tic()
GLM_H <- stan(file = "GLM_Hierarchical.stan", 
            data = data_glm_h,
            chains = 4, 
            iter = 25000, 
            warmup = 5000, 
            thin= 10, 
            seed = 42, 
            init = inits_h,
            control = list(adapt_delta = 0.9,max_treedepth=15),
            algorithm = 'NUTS')
toc()}
save(GLM_H, file="GLM_H.dat")

print(GLM_H, pars=c('beta','sigma2_beta','a','b'))

mcmcModel_GLM_H<-As.mcmc.list(GLM_H)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel_GLM_H[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_GLM_H[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_GLM_H[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_GLM_H[[1]][,4]) # traceplot of beta[4]


## AutoCorrelation Funtion
x11()
acf(mcmcModel_GLM_H[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_GLM_H[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_GLM_H[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_GLM_H[[1]][,4]) #autocorrelation for beta[4]

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_GLM_H)
gelman.diag(mcmcModel_GLM_H, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

graphics.off() #shuts down all open graphics devices

##################### NIMBLE #######################
## Constant values
a0 <- 1/3
b0 <- 1/8

## Initial values
sigma2_0 <- rep(10, p_fix) #real value [4 4 4 4]
tau_0 <- rep(0.1, p_fix)
b_0 <- rep(2,p_fix) #real value [1 2 3 4]
alpha_0 <- c(5,5) #real value [3 8]

###############  HIERARCHICAL MODEL  #####################################
model_Hierarchical<-nimbleCode({
  for(i in 1:N){
    for(j in 1:K){
      a[i,j] <- X[i,j]*beta[j]
    }
    zeta[i] <- sum(a[i,]) # deterministic nodes
    ## likelihood
    Y[i] ~ dpois(exp(zeta[i]))  # stochastic nodes
  }
  
  ## priors
  for(i in 1:K){
    beta[i] ~ dnorm(beta0[i],sd=sigma[i])
    sigmasquared[i] ~  dinvgamma(shape = alpha[1], scale = alpha[2])
    sigma[i] <- pow(sigmasquared[i],0.5)
  }
  alpha[1] ~ dexp(a0)
  alpha[2] ~ dexp(b0)
})

GLMConsts <- list(N=N,K=p_fix)
GLMData <-  list(beta0=beta0,X=X,Y=Y,a0=a0,b0=b0)

GLMInits <- list(beta=b_0,alpha=alpha_0,sigmasquared=sigmasquared)

{tic()
  GLM <- nimbleModel(code = model_Hierarchical,dimensions = list(a=c(N,K)),
                     name = "GLM",constants = GLMConsts, data = GLMData, inits = GLMInits)
  C_GLM <- compileNimble(GLM)
  GLMConf <- configureMCMC(GLM, print = TRUE)
  GLMConf$addMonitors(c("alpha","sigmasquared", "beta"))
  
  GLMMCMC <- buildMCMC(GLMConf)
  C_GLMMCMC <- compileNimble(GLMMCMC, project = GLM)
  {tic()
    Psamples <- runMCMC(C_GLMMCMC, niter=40000, nburnin=20000, thin=10, nchains=4,samplesAsCodaMCMC = TRUE,summary=TRUE)
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


########################################################

################ NON INFORMATIVE PRIOR ################

########################################################
##################################

## Initial values
b_0 <- rep(2,p_fix)

######################## JAGS ###########################

# Compilation Not Hierarchical Model
model_Non_Informative<-
  "model{
  for(i in 1:N){  
    for(j in 1:K){
    a[i,j] <- X[i,j]*beta[j] # deterministic node
    }
    zeta[i] <- sum(a[i,])
    ## likelihood
    Y[i] ~ dpois(exp(zeta[i]))  # stochastic node
  }
  
  ## priors
  for(i in 1:K){
  beta[i] ~ dunif(-100,100) # stochastic node
  } 
}"
model3<-textConnection(model_Non_Informative)

## Compilation phase
{tic()
  JAGS_NI<-jags.model(file=model3,
                   data=list('N'=N, 'X'=X,'K'= p_fix,'Y'=Y),
                   inits=list('beta'=b_0),
                   n.chains=4,
                   n.adapt=5000) #n.adapth=burn-in #(N=10000 -> idem)
  
  ## Sampling phase
  ## Run with CODA
  CODA_NI<-coda.samples(model=JAGS_NI,
                     variable.names=c('beta'), 
                     n.iter=20000,                  #(N=10000 -> idem)
                     thin=10)
  toc()}

## Funtions for manipulating jags.model.objects
list.samplers(JAGS_NI)
coef(JAGS_NI)

## Results by CODA
#head(CODA)
summary(CODA_NI)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA_NI)

## AutoCorrelation Funtion
x11()
acfplot(CODA_NI)


## Effective Sample Size
effectiveSize(CODA_NI)
gelman.diag(CODA_NI)

dev.off()
dev.off()

######################## STAN ###########################



data_glm_ni <-list(N = N, 
                p_fix=p_fix,
                Y=as.vector(Y),
                X=as.matrix(X)
) 


#initialization of the parameters


inits_ni <- function() 
{
  list(
    beta = b_0
  )
}

# run stan model
{tic()
  GLM_NI <- stan(file = "GLM_NonInformative.stan", 
                data = data_glm_ni,
                chains = 4, 
                iter = 25000, 
                warmup = 5000, 
                thin= 10, 
                seed = 42, 
                init = inits_ni,
                algorithm = 'NUTS')
  toc()}
save(GLM_NI, file="GLM_NI.dat")

print(GLM_NI, pars=c('beta'))

mcmcModel_GLM_NI<-As.mcmc.list(GLM_NI)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel_GLM_NI[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_GLM_NI[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_GLM_NI[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_GLM_NI[[1]][,4]) # traceplot of beta[4]


## AutoCorrelation Funtion
x11()
acf(mcmcModel_GLM_NI[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_GLM_NI[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_GLM_NI[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_GLM_NI[[1]][,4]) #autocorrelation for beta[4]

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_GLM_NI)
gelman.diag(mcmcModel_GLM_NI, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

graphics.off() #shuts down all open graphics devices

#################### NIMBLE ################

## Initial values
b_0 <- rep(2,p_fix)

model_Not_Informative<-nimbleCode({
  for(i in 1:N){
    for(j in 1:K){
      a[i,j] <- X[i,j]*beta[j]
    }
    zeta[i] <- sum(a[i,]) # deterministic nodes
    ## likelihood
    Y[i] ~ dpois(exp(zeta[i]))  # stochastic nodes
  }
  
  ## priors
  for(i in 1:K){
    beta[i] ~ dflat()
  }
})

GLMConsts <- list(N=N,K=p_fix)
GLMData <-  list(X=X,Y=Y)

GLMInits <- list(beta=b_0)

{tic()
  GLM <- nimbleModel(code = model_Not_Informative,dimensions = list(a=c(N,K)),
                     name = "GLM",constants = GLMConsts, data = GLMData, inits = GLMInits)
  C_GLM <- compileNimble(GLM)
  GLMConf <- configureMCMC(GLM, print = TRUE)
  GLMConf$addMonitors(c("beta"))
  
  GLMMCMC <- buildMCMC(GLMConf)
  C_GLMMCMC <- compileNimble(GLMMCMC, project = GLM)
  {tic()
    Psamples <- runMCMC(C_GLMMCMC, niter=20000, nburnin=10000, thin=10, nchains=4,samplesAsCodaMCMC = TRUE,summary=TRUE)
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





