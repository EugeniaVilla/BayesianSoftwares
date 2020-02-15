
//

data {
  int<lower=0> N; // number of data
  vector[N] y;    // response
  int<lower=0> p; //number of regressors
  matrix[N, p] X; //matrix of covariates
  
  // Zellner prior parameters
 // real c;  // ONLY THING THAT WORKS IS TO MAKE C A PARAMETER (eith flat prior)
                      //(not as a fixed constant passed as data)--> in this way:
                        // c tands to infinity (meaning that in the posterior mean the prior informations represented by b0 get really low weigth)
                        //--> we give more importance to data and less to prior information
                        
            //otherwise i could fix b0 more close to true values (if we have prior knowledge) 
            // or increase c (c=10000) and use it as a data
  matrix[p, p] B0;
  vector[p] b0; 
  
  // prior parameters of gamma 
  real sigma0;
  real nu0;
}


parameters {
  vector[p] beta;       //regressors paramters (column vector)
  real<lower=0> sigma2;  // standard deviation
  real<lower=0> c; 
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
  beta ~ multi_normal(b0,  sigma2*c*B0); 

  sigma2 ~ inv_gamma(nu0/2, nu0*sigma0/2); // paramters of zellner prior

}

