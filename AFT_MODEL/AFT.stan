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
}

parameters {
  vector[p] beta;                
  real alpha; // Weibull Shape      
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
  beta ~ normal(0, 10);
  alpha ~ exponential(1);
  
  // evaluate likelihood for censored and uncensored subjects
  target += weibull_lpdf(y_o | alpha, lambda_o);
  target += weibull_lccdf(y_m | alpha, lambda_m);
}


// generate posterior quantities of interest
//generated quantities{
//  vector[1000] post_pred_trt;
//  vector[1000] post_pred_pbo;
//  real lambda_trt; 
//  real lambda_pbo; 
//  real hazard_ratio;
//  
//  // generate hazard ratio
//  lambda_trt = exp((beta[1] + beta[2])*alpha ) ;
//  lambda_pbo = exp((beta[1])*alpha ) ;
//  
//  hazard_ratio = exp(beta[2]*alpha ) ;
//  
//  // generate survival times (for plotting survival curves)
//  for(i in 1:1000){
//    post_pred_trt[i] = weibull_rng(alpha,  lambda_trt);
//    post_pred_pbo[i] = weibull_rng(alpha,  lambda_pbo);
//  }
//}
