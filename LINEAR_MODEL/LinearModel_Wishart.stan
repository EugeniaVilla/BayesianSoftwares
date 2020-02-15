
//

data {
  int<lower=0> N; // number of data
  vector[N] y;    // response
  int<lower=0> p; //number of regressors
  matrix[N, p] X; //matrix of covariates
  
  // prior hyperparameters
	matrix[p, p] B0; 
  vector[p] b0; 
  
  // prior parameters of gamma 
  real sigma0;
  real nu0;
  cov_matrix[p] Omega; 
}


parameters {
  vector[p] beta;       //regressors paramters (column vector)
  real<lower=0> sigma2;  // standard deviation
  
  cov_matrix[p] Sigma_beta;
}

transformed parameters 
{
	vector[N] mu;         // mean 
     
	
	for(i in 1:N){
    mu[i] =  X[i,] * beta;
	}


}

model {
  //likelihood:
  y ~ normal(mu, pow(sigma2, 0.5));
  
  //Prior:
 
   Sigma_beta ~ inv_wishart(p+2, Omega); 
   beta ~ multi_normal(b0, Sigma_beta); 
   sigma2 ~ inv_gamma(nu0/2,nu0*sigma0/2 );
}



