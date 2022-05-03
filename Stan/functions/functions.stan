

real poisson_glm_conjprior(vector beta, matrix X, vector m, real lambda) {
    int n = rows(X);
    int p = cols(X);
    vector[n] eta = X * beta;
    return lambda * sum( m .* eta - exp(eta) );
}

real binomial_glm_conjprior(vector beta, matrix X, vector m, real lambda) {
  int n = rows(X);
  vector[n] eta = X * beta;
  return lambda * sum( m .* eta - inv_logit(eta) );
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



real binomial_lognc_irls(vector m, matrix X, real lambda, vector start, real tol, int maxit) {
  real log2pi = 1.83787706641;
  int n       = rows(X);
  int p       = cols(X);
  vector[p] s = start;
  matrix[p,p] XtWX;
  for ( i in 1:maxit ) {
    vector[n]   eta   = X * s;
    vector[n]   mu    = inv_logit(eta);
    matrix[p,n] XtW   = diag_post_multiply(X', mu);
    vector[p]   sold  = s;   // copy s to assess convergence
    
    XtWX = XtW * X;
    s    = mdivide_left_spd( XtWX, XtW * (eta + (m ./ mu - 1) ) );   // update s via WLS with z* = eta + (m - mu) dmu/deta = eta + (m - mu) / mu = eta + (m/mu - 1)
    if ( sqrt( dot_self(s - sold) ) < tol )
      break;
  }
  // return laplace approximation to log nc
  return 0.5 * ( p * log2pi - log(lambda) - log_determinant(XtWX) ) + binomial_glm_conjprior(s, X, m, lambda);
}







