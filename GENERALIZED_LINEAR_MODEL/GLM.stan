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
	vector[p_fix] sigma2_beta;
	vector[p_fix] beta0;
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
	vector[p_fix] beta;        	// regression coefficients 

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
	 	beta[j] ~ normal(beta0[j], pow(sigma2_beta[j], 0.5));

	}
}



