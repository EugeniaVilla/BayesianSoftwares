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

######### DATA GENERATION LINEAR MIXED EFFECT MODEL

########################
# Data Set

## Set the seed for the Random Number Generator
set.seed(431234) # second seed : set.seed(435)

## Create some simulated data

ngr<-3 # number of groups
p_fix<-2  #number of fixed effect parameters (beta)
p_ran<-2 #number of random effect parameters (gamma)

# we suppose that the intercept is modeled separately from the rest (is theta), and is considered as a random effect


# true values of parameters
beta<-seq(1,p_fix,1) #[1 2] # fixed effect regressors
gamma<-matrix( nrow = p_ran, ncol = ngr)
theta<-rep(0,ngr)
for (j in 1:ngr){
  gamma[,j]<-seq(1/j,p_ran/j, 1/j) # random effect regressors
  theta[j]<- j^(1/j)/2 # intercept (a caso)
} 


#	vector[p_fix] sigma2_beta;	        // variances for the prior on beta
# vector[ngr] sigma2_theta;	          // variances for the prior on theta
# matrix[p_ran, ngr] sigma2_gamma;    // (group specific) random effects
sigma2_beta=rep(0,p_fix)
for (i in 1:p_fix)
  sigma2_beta[i]<-1+i*0.6       # variances for the prior on beta

sigma2_theta<-rep(0,ngr)
for (i in 1:ngr)
  sigma2_theta[i]<-1+i^2*0.6 # variances for the prior on theta

sigma2_gamma<-matrix( nrow = p_ran, ncol = ngr) 
for (i in 1:p_ran){
  for (j in 1:ngr){
    sigma2_gamma[i,j]<-j*1.1+i^2*0.6 # (group specific) random effects
  }
}

sigma2_e <- 5


K=750
N_gr=seq(K, K+(ngr-1)*250, 250) # number of data for each group
N<-sum(N_gr)   #total number of data


X<-matrix(nrow=N, ncol=p_fix)
#X[,1]<-rep(1,N) # NO INTERCEPT because the intercept is theta
for(i in 1:p_fix){
  X[,i]<-rnorm(n = N, mean = i^2 ,sd = 2*i)
}

# Matrix X of covariates of random effect
Z<-matrix(nrow=N, ncol=p_ran)
#Z[,1]<-rep(1,N) # NO INTERCEPT
for(i in 1:p_ran){
  Z[,i]<-rnorm(n = N, mean = i^0.5 ,sd = 0.5*i)
}

# G is the matrix that defines the group allocation

G<-matrix(nrow= N, ncol=ngr)


Y_gr=rnorm(n=N_gr[1], mean=X[1:N_gr[1],]%*%beta+Z[1:N_gr[1],]%*%gamma[,1]+rep(theta[1],N_gr[1]),sd=sigma2_e^0.5)
Y<-t(t(Y_gr))

G[1:N_gr[1],]<-rep(0,ngr)
G[1:N_gr[1],1]<-1


ultimo<-N_gr[1]
for (i in 2:ngr){
  inizio=ultimo+1     # first data of the group
  fine=ultimo+N_gr[i] # last data of the group
  Y_gr=rnorm(n=N_gr[i], mean=X[inizio:fine,]%*%beta+Z[inizio:fine,]%*%gamma[,i]+rep(theta[i],N_gr[i]),sd=sigma2_e^0.5)
  Y<-rbind(Y, t(t(Y_gr)))
  
  G[inizio:fine,]<-rep(0, N_gr[i]*ngr)
  G[inizio:fine,i]<-rep(1,N_gr[i])
  ultimo<- fine
}



#########################################################################

####################### NOT HIERARCHICAL PRIOR ##########################

#########################################################################

# initial values
beta0 = rep(0.1, p_fix)
theta0=rep(0.1, ngr)
gamma0=matrix(0.1, p_ran, ngr)

sigma2_beta0=rep(10, p_fix)
sigma2_theta0=rep(10,ngr)
sigma2_gamma0=matrix(10, p_ran,ngr)
sigma2_e0=10

################# JAGS #########################

# Compilation Model phase
model_simplified<-
  "model{
  #Likelihood     
	for (i in 1:N){
  for(s in 1:ngr){
    a[i,s] <- Z[i,1:p_ran]%*%t(t(gamma[1:p_ran,s]))*G[i,s] + theta[s]*G[i,s]
  }
  for(k in 1:p_fix){
      c[i,k] <- X[i,k]*beta[k]
    }
    mu[i]<- sum(a[i,])+sum(c[i,])
		Y[i,1] ~ dnorm(mu[i],tau)
	} 
  # Priors
  # beta
	for (j in 1:p_fix){
		beta[j] ~ dnorm(0.0, tau_beta[j])
		sigma2_beta[j] <- pow(tau_beta[j], -1)
		tau_beta[j] ~ dgamma(2., 10.)
	 }
  # theta
	for (j in 1:ngr){
	 	theta[j] ~ dnorm(0.0, tau_theta[j])
	 	sigma2_theta[j] <- pow(tau_theta[j], -1)
		tau_theta[j] ~ dgamma(2., 10.)
	}
	# gamma
	for(j in 1:p_ran){
	for (i in 1:ngr){
	 	gamma[j,i] ~ dnorm(0.0, tau_gamma[j,i])
	 	sigma2_gamma[j,i] <- pow(tau_gamma[j,i], -1)
	  tau_gamma[j,i] ~ dgamma(2., 10.)
	}
	}
	sigma2_e <- pow(tau,-1)
	tau ~ dgamma(2., 10.)
}"
model1<-textConnection(model_simplified)

## Compilation phase
{tic()
  JAGS1<-jags.model(file=model1,
                  data=list('N'=N,'Y'=Y,'ngr'= ngr,'p_ran'=p_ran,'Z'=Z,'G'=G,'p_fix'=p_fix,'X'=X),
                  inits=list('gamma'=gamma0,'beta'=beta0,'theta'=theta0,'tau'=1/sigma2_e0,'tau_beta'=1/sigma2_beta0,'tau_theta'=1/sigma2_theta0,'tau_gamma'=1/sigma2_gamma0),
                  n.chains=4,
                  n.adapt=5000) #n.adapth=burn-in
## Sampling phase
## Run with CODA
CODA1<-coda.samples(model=JAGS1,
                    variable.names=c('gamma','beta','theta'), 
                    n.iter=5000,
                    thin=10) 
  toc()
}

## Funtions for manipulating jags.model.objects
list.samplers(JAGS1)
coef(JAGS1)

## Results by CODA
#head(CODA1)
summary(CODA1)

beta
gamma
theta

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA1)

## AutoCorrelation Funtion
x11()
acfplot(CODA1)


## Effective Sample Size
effectiveSize(CODA1)
gelman.diag(CODA1)


dev.off()
dev.off()

################# STAN #########################

data_lmm <-list(N = N, 
                p_fix=p_fix,
                ngr=ngr,
                p_ran=p_ran,
                Y=as.vector(Y),
                X=as.matrix(X),
                G=as.matrix(G),
                Z=as.matrix(Z)     # <3
) 


#initialization of the parameters

inits <- function() 
{
  list(
    beta =beta0, 
    theta=theta0,
    gamma=gamma0, 
    
    sigma2_beta=sigma2_beta0,
    sigma2_theta=sigma2_theta0,
    sigma2_gamma=sigma2_gamma0,
    sigma2_e=sigma2_e0)
}

# run stan model
{tic()
LMM <- stan(file = "LinearMixedModel.stan", 
            data = data_lmm,
            chains = 1, 
            iter = 15000, 
            warmup = 5000, 
            thin= 5, 
            seed = 42, 
            init = inits,
            algorithm = 'NUTS')
toc()}

save(LMM, file="LMM.dat")
load("LMM.dat")

print(LMM, pars=c('beta','gamma','theta','sigma2_e'))

mcmcModel_LMM<-As.mcmc.list(LMM)

## Diagnostic by CODA
## Traceplot and Density plot
#order in which parameters appears in mcmcModel:
#betas
#thetas
#gamma[1,i]
#gamma[2,i]
# sigmas in the same order

x11()
coda::traceplot(mcmcModel_LMM[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_LMM[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_LMM[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_LMM[[1]][,4]) # traceplot of beta[4]


## AutoCorrelation Funtion
x11()
acf(mcmcModel_LMM[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_LMM[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_LMM[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_LMM[[1]][,4]) #autocorrelation for beta[4]

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_LMM)
gelman.diag(mcmcModel_LMM, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

graphics.off() #shuts down all open graphics devices

############### NIMBLE ######################################

# initial values
beta0 = rep(0.1, p_fix)
theta0=rep(0.1, ngr)
gamma0=matrix(0.1, p_ran, ngr)

sigma2_beta0=rep(10, p_fix)
sigma2_theta0=rep(10,ngr)
sigma2_gamma0=matrix(10, p_ran,ngr)
sigma2_e0=10

linearMixed_NH_Code <- nimbleCode({
  #Likelihood     
  for (i in 1:N){
    for(s in 1:ngr){
      a[i,s] <- (Z[i,1:p_ran]%*%t(t(gamma[1:p_ran,s]))*G[i,s] + theta[s]*G[i,s])[1,1]
    }
    for(k in 1:p_fix){
      c[i,k] <- X[i,k]*beta[k]
    }
    mu[i]<- sum(a[i,])+sum(c[i,])
    Y[i,1] ~ dnorm(mu[i],sd = sigma_e)
  } 
  # Priors
  # beta
  for (j in 1:p_fix){
    beta[j] ~ dnorm(0.0, sd=sigma_beta[j])
    sigma2_beta[j] ~ dinvgamma(2., 10.)
    sigma_beta[j] <- pow(sigma2_beta[j],0.5)
  }
  # theta
  for (j in 1:ngr){
    theta[j] ~ dnorm(0.0, sd = sigma_theta[j])
    sigma2_theta[j] ~ dinvgamma(2., 10.)
    sigma_theta[j] <- pow(sigma2_theta[j], 0.5)
  }
  # gamma
  for(j in 1:p_ran){
    for (i in 1:ngr){
      gamma[j,i] ~ dnorm(0.0, sd = sigma_gamma[j,i])
      sigma2_gamma[j,i] ~ dinvgamma(2., 10.)
      sigma_gamma[j,i]  <- pow(sigma2_gamma[j,i],0.5)
    }
  }
  sigma2_e ~ dinvgamma(2, 10.)
  sigma_e <- pow(sigma2_e,0.5)
})


linearmixedConsts <- list(N=N,ngr=ngr,p_ran=p_ran,p_fix=p_fix)
linearmixedInits <-  list(gamma=gamma0,beta=beta0,theta=theta0,sigma2_e=sigma2_e0,
                          sigma2_beta=sigma2_beta0,sigma2_theta=sigma2_theta0,sigma2_gamma=sigma2_gamma0)
linearmixedData <- list(X=X,Y=Y,G=G,Z=Z)

{tic()
  linearmixedNH <- nimbleModel(code = linearMixed_NH_Code,dimensions = list(a = c(N, ngr),c=c(N,p_fix)),
                               name = "linearmixedNH",constants = linearmixedConsts, data = linearmixedData, inits = linearmixedInits)
  Clinearmixed <- compileNimble(linearmixedNH)
  linearmixedConf <- configureMCMC(linearmixedNH, print = TRUE)
  linearmixedConf$addMonitors(c("gamma", "beta","theta"))
  linearmixedMCMC <- buildMCMC(linearmixedConf)
  ClinearmixedMCMC <- compileNimble(linearmixedMCMC, project = linearmixedNH)
  {tic()
    Psamples <- runMCMC(ClinearmixedMCMC, niter=10000, nburnin=5000, thin=10,
                        nchains=2,samplesAsCodaMCMC = TRUE,summary=TRUE)
    toc()}
  toc()}

Psamples$summary

gelman.diag(Psamples[[1]], confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
effectiveSize(Psamples[[1]])


#########################################################################

########################### HIERARCHICAL PRIOR ##########################

#########################################################################

# initial values
beta0 = rep(0.1, p_fix)
theta0=rep(0.1, ngr)
gamma0=matrix(0.1, p_ran, ngr)

sigma2_beta0=rep(10, p_fix)
sigma2_theta0=rep(10,ngr)
sigma2_gamma0=matrix(10, p_ran,ngr)
sigma2_e0=10


#initial values for additional parameters
a0=1
b0=1

ab0=1
bb0=1

at0=1
bt0=1

ag0=1
bg0=1

# fixed values of hyperparamters
a_hat=1/2
b_hat=1/10

################# JAGS #########################

# Compilation Model phase
model_hier<-
  "model{
  #Likelihood     
	for (i in 1:N){
  for(s in 1:ngr){
    aa[i,s] <- Z[i,1:p_ran]%*%t(t(gamma[1:p_ran,s]))*G[i,s] + theta[s]*G[i,s]
  }
  for(k in 1:p_fix){
      c[i,k] <- X[i,k]*beta[k]
    }
    mu[i]<- sum(aa[i,])+sum(c[i,])
		Y[i,1] ~ dnorm(mu[i],tau)
	} 
  # Priors
  # beta
	for (j in 1:p_fix){
		beta[j] ~ dnorm(0.0, tau_beta[j])
		sigma2_beta[j] <- pow(tau_beta[j], -1)
		tau_beta[j] ~ dgamma(ab, bb)
	 }
  # theta
	for (j in 1:ngr){
	 	theta[j] ~ dnorm(0.0, tau_theta[j])
	 	sigma2_theta[j] <- pow(tau_theta[j], -1)
		tau_theta[j] ~ dgamma(at, bt)
	}
	# gamma
	for(j in 1:p_ran){
	for (i in 1:ngr){
	 	gamma[j,i] ~ dnorm(0.0, tau_gamma[j,i])
	 	sigma2_gamma[j,i] <- pow(tau_gamma[j,i], -1)
	  tau_gamma[j,i] ~ dgamma(ag, bg)
	}
	}
	sigma2_e <- pow(tau,-1)
	tau ~ dgamma(a, b)
	a ~ dexp(a_hat)
	b ~ dexp(b_hat)
	ag ~ dexp(a_hat)
	bg ~ dexp(b_hat)
	at ~ dexp(a_hat)
	bt ~ dexp(b_hat)
	ab ~ dexp(a_hat)
	bb ~ dexp(b_hat)
}"
model2<-textConnection(model_hier)

## Compilation phase
{tic()
  JAGS2<-jags.model(file=model2,
                    data=list('N'=N,'Y'=Y,'ngr'= ngr,'p_ran'=p_ran,'Z'=Z,'G'=G,'p_fix'=p_fix,'X'=X,'a_hat'=a_hat,'b_hat'=b_hat),
                    inits=list('gamma'=gamma0,'beta'=beta0,'theta'=theta0,'tau'=1/sigma2_e0,'tau_beta'=1/sigma2_beta0,'tau_theta'=1/sigma2_theta0,'tau_gamma'=1/sigma2_gamma0,
                               'a'=a0,'b'=b0,'ab'=ab0,'bb'=bb0,'ag'=ag0,'bg'=bg0,'at'=at0,'bt'=bt0),
                    n.chains=4,
                    n.adapt=5000) #n.adapth=burn-in
  ## Sampling phase
  ## Run with CODA
  CODA2<-coda.samples(model=JAGS2,
                      variable.names=c('gamma','beta','theta','sigma2_e'), 
                      n.iter=5000,
                      thin=10) 
  toc()
}
## Funtions for manipulating jags.model.objects
list.samplers(JAGS2)
coef(JAGS2)

## Results by CODA
#head(CODA1)
summary(CODA2)

beta
gamma
theta

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA2)

## AutoCorrelation Funtion
x11()
acfplot(CODA2)


## Effective Sample Size
effectiveSize(CODA2)
gelman.diag(CODA2)

dev.off()
dev.off()

################# STAN #########################

data_lmm_H <-list(N = N, 
                p_fix=p_fix,
                ngr=ngr,
                p_ran=p_ran,
                Y=as.vector(Y),
                X=as.matrix(X),
                G=as.matrix(G),
                Z=as.matrix(Z),     
                a_hat=a_hat,
                b_hat=b_hat
) 


#initialization of the parameters

inits_H <- function() 
{
  list(
    beta =beta0, 
    theta=theta0,
    gamma=gamma0, 
    
    sigma2_beta=sigma2_beta0,
    sigma2_theta=sigma2_theta0,
    sigma2_gamma=sigma2_gamma0,
    sigma2_e=sigma2_e0,
    
    a=a0,
    b=b0,
    ab=ab0,
    bb=bb0,
    at=at0,
    bt=bt0,
    ag=ag0,
    bg=bg0
    )
}

# run stan model
{tic()
LMM_H <- stan(file = "LMM_Hierarchical.stan", 
            data = data_lmm_H,
            chains = 1, 
            iter = 10000, 
            warmup = 3000, 
            thin= 5, 
            seed = 42, 
            init = inits_H,
            algorithm = 'NUTS')
toc()}

save(LMM_H, file="LMM.dat")
load("LMM_H.dat")

print(LMM_H, pars=c('beta','gamma','theta','sigma2_e'))

x11()
rstan::traceplot(LMM_H, pars='beta') 

x11()
rstan::traceplot(LMM_H, pars='theta') 

x11()
rstan::traceplot(LMM_H, pars='gamma') 



mcmcModel_LMM_H<-As.mcmc.list(LMM_H)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel_LMM_H[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_LMM_H[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_LMM_H[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_LMM_H[[1]][,4]) # traceplot of beta[4]


## AutoCorrelation Funtion
x11()
acf(mcmcModel_LMM_H[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_LMM_H[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_LMM_H[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_LMM_H[[1]][,4]) #autocorrelation for beta[4]

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_LMM_H)
gelman.diag(mcmcModel_LMM_H, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

graphics.off() #shuts down all open graphics devices

############### NIMBLE ######################################

# initial values
beta0 = rep(0.1, p_fix)
theta0=rep(0.1, ngr)
gamma0=matrix(0.1, p_ran, ngr)

sigma2_beta0=rep(10, p_fix)
sigma2_theta0=rep(10,ngr)
sigma2_gamma0=matrix(10, p_ran,ngr)
sigma2_e0=10


#initial values for additional parameters
a0=1
b0=1

ab0=1
bb0=1

at0=1
bt0=1

ag0=1
bg0=1

# fixed values of hyperparamters
a_hat=1/2
b_hat=1/10
#modello
linearMixed_H_Code <- nimbleCode({
  #Likelihood     
  for (i in 1:N){
    for(s in 1:ngr){
      a[i,s] <- (Z[i,1:p_ran]%*%t(t(gamma[1:p_ran,s]))*G[i,s] + theta[s]*G[i,s])[1,1]
    }
    for(k in 1:p_fix){
      c[i,k] <- X[i,k]*beta[k]
    }
    mu[i]<- sum(a[i,])+sum(c[i,])
    Y[i,1] ~ dnorm(mu[i],sd = sigma_e)
  } 
  # Priors
  # beta
  for (j in 1:p_fix){
    beta[j] ~ dnorm(0.0, sd=sigma_beta[j])
    sigma2_beta[j] ~ dinvgamma(ab,bb)
    sigma_beta[j] <- pow(sigma2_beta[j],0.5)
  }
  # theta
  for (j in 1:ngr){
    theta[j] ~ dnorm(0.0, sd = sigma_theta[j])
    sigma2_theta[j] ~ dinvgamma(at,bt)
    sigma_theta[j] <- pow(sigma2_theta[j], 0.5)
  }
  # gamma
  for(j in 1:p_ran){
    for (i in 1:ngr){
      gamma[j,i] ~ dnorm(0.0, sd = sigma_gamma[j,i])
      sigma2_gamma[j,i] ~ dinvgamma(ag,bg)
      sigma_gamma[j,i]  <- pow(sigma2_gamma[j,i],0.5)
    }
  }
  sigma2_e ~ dinvgamma(a,b)
  sigma_e <- pow(sigma2_e,0.5)
  a ~ exp(a_hat)
  b ~ exp(b_hat)
  ag ~ exp(a_hat)
  bg ~ exp(b_hat)
  at ~ exp(a_hat)
  bt ~ exp(b_hat)
  ab ~ exp(a_hat)
  bb ~ exp(b_hat)
})


linearmixedConsts <- list(N=N,ngr=ngr,p_ran=p_ran,p_fix=p_fix)
linearmixedInits <-  list(gamma=gamma0,beta=beta0,theta=theta0,sigma2_e=sigma2_e0,
                          sigma2_beta=sigma2_beta0,sigma2_theta=sigma2_theta0,sigma2_gamma=sigma2_gamma0,
                          a=a0,b=b0,ab=ab0,bb=bb0,ag=ag0,bg=bg0,at=at0,bt=bt0)
linearmixedData <- list(X=X,Y=Y,G=G,Z=Z)

{tic()
  linearmixedH <- nimbleModel(code = linearMixed_H_Code,dimensions = list(a = c(N, ngr),c=c(N,p_fix)),
                              name = "linearmixedH",constants = linearmixedConsts, data = linearmixedData, inits = linearmixedInits)
  Clinearmixed <- compileNimble(linearmixedH)
  linearmixedConf <- configureMCMC(linearmixedH, print = TRUE)
  linearmixedConf$addMonitors(c("gamma", "beta","theta"))
  linearmixedMCMC <- buildMCMC(linearmixedConf)
  ClinearmixedMCMC <- compileNimble(linearmixedMCMC, project = linearmixedH)
  {tic()
    Psamples <- runMCMC(ClinearmixedMCMC, niter=100000, nburnin=50000, thin=10,
                        nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
    toc()}
  toc()}

Psamples$summary

gelman.diag(Psamples[[1]], confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
effectiveSize(Psamples[[1]])

########## Traceplot
x11()
traceplot(Psamples[["samples"]][,7], xlab = "iteration", ylab = expression(alpha),main="sigma2_beta[1] traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,5],xlab = "iteration", ylab = expression(sbeta),main="sigmasq traceplot")
x11()
traceplot(Psamples[["samples"]][["chain1"]][,6],xlab = "iteration", ylab = expression(sbeta),main="tau traceplot")

x11()   #-----------   AUTOCORRELATION
acf(Psamples[["samples"]][,7],main="sigma2_beta[1] autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,5],main="sigmasq autocorrelation")
x11()
acf(Psamples[["samples"]][["chain1"]][,6],main="tau autocorrelation")

#############################################################

################### NON INFORMATIVE PRIOR ###################

#############################################################


# initial values
beta0 = rep(0.1, p_fix)
theta0=rep(0.1, ngr)
gamma0=matrix(0.1, p_ran, ngr)

sigma2_beta0=rep(10, p_fix)
sigma2_theta0=rep(10,ngr)
sigma2_gamma0=matrix(10, p_ran,ngr)
sigma2_e0=10

################# JAGS #########################

# Compilation Model phase
model_not<-
  "model{
  #Likelihood     
	for (i in 1:N){
  for(s in 1:ngr){
    a[i,s] <- Z[i,1:p_ran]%*%t(t(gamma[1:p_ran,s]))*G[i,s] + theta[s]*G[i,s]
  }
  for(k in 1:p_fix){
      c[i,k] <- X[i,k]*beta[k]
    }
    mu[i]<- sum(a[i,])+sum(c[i,])
		Y[i,1] ~ dnorm(mu[i],tau)
	} 
  # Priors
  # beta
	for (j in 1:p_fix){
		beta[j] ~ dunif(-100, 100)
	 }
  # theta
	for (j in 1:ngr){
	 	theta[j] ~ dunif(-100, 100)
	}
	# gamma
	for(j in 1:p_ran){
	for (i in 1:ngr){
	 	gamma[j,i] ~ dunif(-100, 100)
	}
	}
	sigma2_e <- pow(tau,-1)
	tau ~ dunif(0.0001, 10000)
}"
model3<-textConnection(model_not)

## Compilation phase
{tic()
  JAGS3<-jags.model(file=model3,
                    data=list('N'=N,'Y'=Y,'ngr'= ngr,'p_ran'=p_ran,'Z'=Z,'G'=G,'p_fix'=p_fix,'X'=X),
                    inits=list('gamma'=gamma0,'beta'=beta0,'theta'=theta0,'tau'=1/sigma2_e0),
                    n.chains=4,
                    n.adapt=5000) #n.adapth=burn-in
  ## Sampling phase
  ## Run with CODA
  CODA3<-coda.samples(model=JAGS3,
                      variable.names=c('gamma','beta','theta','sigma2_e'), 
                      n.iter=5000,
                      thin=10) 
  toc()
}
## Funtions for manipulating jags.model.objects
list.samplers(JAGS3)
coef(JAGS3)

## Results by CODA
#head(CODA3)
summary(CODA3)

beta
gamma
theta

## Diagnostic by CODA
## Traceplot and Density plot
x11()
plot(CODA3)

## AutoCorrelation Funtion
x11()
acfplot(CODA3)


## Effective Sample Size
effectiveSize(CODA3)
gelman.diag(CODA3)

dev.off()
dev.off()

################# STAN #########################

data_lmm_NI <-list(N = N, 
                p_fix=p_fix,
                ngr=ngr,
                p_ran=p_ran,
                Y=as.vector(Y),
                X=as.matrix(X),
                G=as.matrix(G),
                Z=as.matrix(Z)     
) 


#initialization of the parameters

inits_NI <- function() 
{
  list(
    beta =beta0, 
    theta=theta0,
    gamma=gamma0, 
    
    sigma2_e=sigma2_e0)
}

# run stan model
{tic()
LMM_NI <- stan(file = "LMM_Non_Informative.stan", 
            data = data_lmm_NI,
            chains = 1, 
            iter = 10000, 
            warmup = 3000, 
            thin= 5, 
            seed = 42, 
            init = inits_NI,
            algorithm = 'NUTS')
toc()}

save(LMM_NI, file="LMM_NI.dat")
load("LMM_NI.dat")

print(LMM_NI, pars=c('beta','gamma','theta','sigma2_e'))

mcmcModel_LMM_NI<-As.mcmc.list(LMM_NI)

## Diagnostic by CODA
## Traceplot and Density plot
x11()
coda::traceplot(mcmcModel_LMM_NI[[1]][,1]) # traceplot of beta[1]
x11()
coda::traceplot(mcmcModel_LMM_NI[[1]][,2]) # traceplot of beta[2]
x11()
coda::traceplot(mcmcModel_LMM_NI[[1]][,3]) # traceplot of beta[3]
x11()
coda::traceplot(mcmcModel_LMM_NI[[1]][,4]) # traceplot of beta[4]


## AutoCorrelation Funtion
x11()
acf(mcmcModel_LMM_NI[[1]][,1]) #autocorrelation for beta[1]
x11()
acf(mcmcModel_LMM_NI[[1]][,2]) #autocorrelation for beta[2]
x11()
acf(mcmcModel_LMM_NI[[1]][,3]) #autocorrelation for beta[3]
x11()
acf(mcmcModel_LMM_NI[[1]][,4]) #autocorrelation for beta[4]

## Effective Sample Size and Rhat
effectiveSize(mcmcModel_LMM_NI)
gelman.diag(mcmcModel_LMM_NI, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

graphics.off() #shuts down all open graphics devices

####################### NIMBLE ###########################


# initial values
beta0 = rep(0.1, p_fix)
theta0=rep(0.1, ngr)
gamma0=matrix(0.1, p_ran, ngr)

sigma2_beta0=rep(10, p_fix)
sigma2_theta0=rep(10,ngr)
sigma2_gamma0=matrix(10, p_ran,ngr)
sigma2_e0=10

linearMixed_NI_Code <- nimbleCode({
  #Likelihood     
  for (i in 1:N){
    for(s in 1:ngr){
      a[i,s] <- (Z[i,1:p_ran]%*%t(t(gamma[1:p_ran,s]))*G[i,s] + theta[s]*G[i,s])[1,1]
    }
    for(k in 1:p_fix){
      c[i,k] <- X[i,k]*beta[k]
    }
    mu[i]<- sum(a[i,])+sum(c[i,])
    Y[i,1] ~ dnorm(mu[i],sd = sigma_e)
  } 
  # Priors
  # beta
  for (j in 1:p_fix){
    beta[j] ~ dunif(-100,100)
  }
  # theta
  for (j in 1:ngr){
    theta[j] ~ dunif(-100,100)
  }
  # gamma
  for(j in 1:p_ran){
    for (i in 1:ngr){
      gamma[j,i] ~ dunif(-100,100)
    }
  }
  sigma2_e ~ dunif(0,1000)
  sigma_e <- pow(sigma2_e,0.5)
})


linearmixedConsts <- list(N=N,ngr=ngr,p_ran=p_ran,p_fix=p_fix)
linearmixedInits <-  list(gamma=gamma0,beta=beta0,theta=theta0)
linearmixedData <- list(X=X,Y=Y,G=G,Z=Z)

{tic()
  linearmixedNI <- nimbleModel(code = linearMixed_NI_Code,dimensions = list(a = c(N, ngr),c=c(N,p_fix)),
                               name = "linearmixedNI",constants = linearmixedConsts, data = linearmixedData, inits = linearmixedInits)
  Clinearmixed <- compileNimble(linearmixedNI)
  linearmixedConf <- configureMCMC(linearmixedNI, print = TRUE)
  linearmixedConf$addMonitors(c("gamma", "beta","theta"))
  linearmixedMCMC <- buildMCMC(linearmixedConf)
  ClinearmixedMCMC <- compileNimble(linearmixedMCMC, project = linearmixedNI)
  {tic()
    Psamples <- runMCMC(ClinearmixedMCMC, niter=1000, nburnin=500, thin=10,
                        nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
    toc()}
  toc()}

Psamples$summary

gelman.diag(Psamples[[1]], confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)
effectiveSize(Psamples[[1]])

