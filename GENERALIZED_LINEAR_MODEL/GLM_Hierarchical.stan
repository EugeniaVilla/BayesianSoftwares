//#######################################################
//### GLM - POISSON REGRESSION 	      ###################
//#######################################################

///////////////////////// DATA /////////////////////////////////////
data {
	int<lower = 0> N;       // number of data
	int<lower = 0> p_fix;   // number of covariates, fixed effects
	
	int<lower = 0> Y[N];  	// response vector
	matrix[N, p_fix] X;   	// design matrix (fixed effects)
	
	//hyperparamters
	vector[p_fix] beta0;
	real<lower=0> a0;
	real<lower=0> b0;
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
	vector[p_fix] beta;        	// regression coefficients 
	vector[p_fix] sigma2_beta;	// variances for the prior on beta
	
	real<lower=0> a;
	real<lower=0> b;
}

//////////////////// TRANSFORMED PARAMETERS /////////////////////////////////
transformed parameters 
{
	vector[N] mu;
	for(i in 1:N){
	  mu[i] = exp(row(X, i) * beta);
	}
}

////////////////// MODEL ////////////////////////
model {

	// Likelihood     
	for (s in 1:N)
	{
		Y[s] ~ poisson(mu[s]);  
	} 

	for (j in 1:p_fix) 
	{
	 	beta[j] ~ normal(0.0, pow(sigma2_beta[j], 0.5));
		sigma2_beta[j] ~ inv_gamma(a, b);
	}
	
  a ~ exponential(a0);
  b ~ exponential(b0);
}

