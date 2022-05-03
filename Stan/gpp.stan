

data {
  int<lower=0>          n;             // number of observations in current data
  int<lower=0>          p;             // number of covariates
  int<lower=0>          y[n];          // integer vector giving responses of current data
  matrix[n, p]          X;             // design matrix for current data
  real<lower=0,upper=1> lambda;        // power prior parameter
  vector[p]             betahat;       // p-dimensional prior mean
  matrix[p,p]           cov_betahat;   // pxp prior covariance
  int<lower=0,upper=1>  post_sample;   // whether to sample from prior or posterior
}
// create power prior cov mtx
transformed data{
  matrix[p,p] cov_gpp = cov_betahat / lambda;
}
// single parameter: regression coefficient
parameters {
  vector[p] beta;
}

// sample from prior if n == 0; else sample from posterior
//    beta ~ N(betahat, cov_gpp) a priori
model {
  beta ~ multi_normal(betahat, cov_gpp);
  if ( post_sample == 1 ) 
    target += poisson_log_glm_lpmf(y | X, 0, beta); 
}

