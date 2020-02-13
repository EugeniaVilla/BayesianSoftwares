data {
  int p; // number of beta parameters
  
  // data for censored subjects
  int N_m;
  matrix[N_m,p] X_m;
  vector[N_m] y_m;
  
  // data for observed subjects
  int N_o;
  matrix[N_o,p] X_o;
  real y_o[N_o];
  
  // hyperparameters
  real mu_0;
  real<lower=0> nu_0;
  real<lower=0> sigma2_0;
  real<lower=0> c;
  real<lower=0> d;
}

parameters {
  vector[p] beta;                
  real alpha; // Weibull Shape  
  real mu;
  real<lower=0> sigma2;
}

transformed parameters{
  // model Weibull rate as function of covariates
  vector[N_m] lambda_m;
  vector[N_o] lambda_o;
  
  // standard weibull AFT re-parameterization
  lambda_m = exp((X_m*beta))/pow(log(2),(1/alpha));
  lambda_o = exp((X_o*beta))/pow(log(2),(1/alpha));
}

model {
  beta ~ normal(mu, pow(sigma2,0.5));
  alpha ~ gamma(c,d);
  
  
  mu ~ normal(mu_0,pow(sigma2,0.5));
  sigma2~inv_gamma(nu_0/2,nu_0*sigma2_0/2);

  
  
  // evaluate likelihood for censored and uncensored subjects
  target += weibull_lpdf(y_o | alpha, lambda_o);
  target += weibull_lccdf(y_m | alpha, lambda_m);
}



