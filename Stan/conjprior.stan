

data {
  int<lower=0>           n;           // current data sample size
  int<lower=0>           p;           // number of covariates (incl. intercept)
  int<lower=0>           y[n];        // vector of integers for poisson responses in current data
  matrix[n, p]           X;           // n x p design matrix for current data
  real<lower=0,upper=1>  lambda;      // precision parameter in conjugate prior
  real<lower=0>          mu0[n];      // prior guess for responses in conjugate prior
  int<lower=0,upper=1>   post_sample; // 1 samples from posterior; 0 samples from prior
}

transformed data{
  vector[n] mu0vec = to_vector(mu0);
}

// single parameter: regression coefficient
parameters {
  vector[p] beta;
}

// sample from power prior if n == 0; else sample from posterior
model {
  vector[n] eta = X * beta;
  target += lambda * ( sum(mu0vec .* eta - exp(eta) ) );
  if ( post_sample == 1 )
    target += poisson_log_glm_lpmf(y | X, 0, beta);
}

