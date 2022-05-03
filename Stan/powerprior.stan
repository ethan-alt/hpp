

data {
  int<lower=0>          n0;          // historical data sample size
  int<lower=0>          n;           // current data sample size
  int<lower=0>          p;           // number of covariates (incl. intercept)
  int<lower=0>          y0[n0];      // vector of integers for poisson responses in historical data
  int<lower=0>          y[n];        // vector of integers for poisson responses in current data
  matrix[n0, p]         X0;          // n0 x p design matrix for historical data
  matrix[n, p]          X;           // n x p design matrix for current data
  real<lower=0,upper=1> lambda;      // power prior parameter
  int<lower=0,upper=1>  post_sample; // indicator for whether to sample from posterior (1 = posterior; 0 = prior)
}

// single parameter: regression coefficient
parameters {
  vector[p] beta;
}

// sample from power prior if n == 0; else sample from posterior
model {
  target += lambda * poisson_log_glm_lpmf(y0 | X0, 0, beta);
  if ( post_sample == 1 )
    target += poisson_log_glm_lpmf(y | X, 0, beta);
}

