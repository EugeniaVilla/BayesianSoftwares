
//

data {
  int<lower=0> N; // number of data
  vector[N] y;    // response
  int<lower=0> p; //number of regressors
  matrix[N, p] X; //matrix of covariates
  

}


parameters {
  vector[p] beta;       //regressors paramters (column vector)

  real<lower=0> sigma2;
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

}
