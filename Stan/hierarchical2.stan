functions {
  real poisson_glm_conjprior(vector beta, matrix X, vector m, real lambda) {
    int n = rows(X);
    int p = cols(X);
    vector[n] eta = X * beta;
    return lambda * sum( m .* eta - exp(eta) );
  }
  
  real poisson_lognc_irls(vector m, matrix X, real lambda, vector start, real tol, int maxit) {
    real log2pi = 1.83787706641;
    int n       = rows(X);
    int p       = cols(X);
    vector[p] s = start;
    matrix[p,p] XtWX;
    for ( i in 1:maxit ) {
      vector[n]   eta   = X * s;
      vector[n]   mu    = exp(eta);
      matrix[p,n] XtW   = diag_post_multiply(X', mu);
      vector[p]   sold  = s;   // copy s to assess convergence
      
      XtWX = XtW * X;
      s    = mdivide_left_spd( XtWX, XtW * (eta + (m ./ mu - 1) ) );   // update s via WLS with z* = eta + (m - mu) dmu/deta = eta + (m - mu) / mu = eta + (m/mu - 1)
      if ( sqrt( dot_self(s - sold) ) < tol )
        break;
    }
    // return laplace approximation to log nc
    return 0.5 * ( p * log2pi - log(lambda) - log_determinant(XtWX) ) + poisson_glm_conjprior(s, X, m, lambda);
  }
}

data {
  int<lower=0>           n;           // current data sample size
  int<lower=0>           p;           // number of covariates (incl. intercept)
  int<lower=0>           y[n];        // vector of integers for poisson responses in current data
  matrix[n, p]           X;           // n x p design matrix for current data
  real<lower=0,upper=1>  lambda;      // precision parameter in conjugate prior
  real<lower=0>          mu0[n];      // prior guess for responses in conjugate prior
  real<lower=0>          lambda0[n];  // vector of precision parameters for hierarchical model
  vector[p]              start;       // starting value for IRLS
  int<lower=0>           maxit;       // max number of iterations for IRLS
  real<lower=0>          tol;         // tolerance for convergence for IRLS
  int<lower=0,upper=1>   post_sample; // 1 samples from posterior; 0 samples from prior
}

transformed data{
  real shape[n]   = to_array_1d( to_vector(lambda0) .* to_vector(mu0) );
  vector[n] yvec  = to_vector(y);
  real lambdapost = 1 + lambda;
}

// single parameter: regression coefficient
parameters {
  vector[p] beta;
  vector[n] m;
}

// sample from power prior if n == 0; else sample from posterior
model {
  m ~ gamma(shape, lambda0);                                       // prior for m is gamma
  target += -poisson_lognc_irls(m, X, lambda, start, tol, maxit);  // subtract log nc of conj prior via laplace
  
  // if post_sample == 0, sample from prior
  if ( post_sample == 0 )
    target += poisson_glm_conjprior(beta, X, m, lambda);
  else  // otherwise, sample from posterior with m* = (lambda * m + y) / (1 + lambda)
    target += poisson_glm_conjprior(beta, X, (lambda * m + yvec) / lambdapost, lambdapost);
}

