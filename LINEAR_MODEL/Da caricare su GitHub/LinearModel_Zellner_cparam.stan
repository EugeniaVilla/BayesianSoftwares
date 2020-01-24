
//

data {
  int<lower=0> N; // number of data
  vector[N] y;    // response
  int<lower=0> p; //number of regressors
  matrix[N, p] X; //matrix of covariates
  
  // Zellner prior parameters
 // real c; //UNICA COSA CHE FUNZIONA é RENDERE QUESTO UN PARAMETRO (con flat prior)  
                      //(non una costante prefissata passata come dato) --> allora così:
            //c tende ad inifinito (vuol dire dare pochissimo peso a b0 nella posterior mean) e sigma2->6
            //diamo più importanza ai dati e meno alle prior information
            
            // oppure mettere b0 più simili alle beta finali 
            // oppure aumentare c (c=10000) e lasciarlo come dato da comunque buoni risultati 
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

