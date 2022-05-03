functions{
  real glm_lp(vector y, vector beta, real phi, matrix X, int dist, int link, vector offset) {
    int n    = rows(y);
    real tau = inv(phi);
    vector[n] eta = X * beta + offset;
    vector[n] mu;
    if      (link == 1) mu = eta;                      // identity link
    else if (link == 2) mu = exp(eta);                 // log link
    else if (link == 3) mu = inv_logit(eta);           // logit link
    else if (link == 4) mu = inv(eta);                 // inverse link
    else if (link == 5) mu = Phi_approx(eta);          // probit link
    else if (link == 6) mu = atan(eta) / pi() + 0.5;   // cauchit link
    else if (link == 7) mu = inv_cloglog(eta);         // complementary log-log link
    else if (link == 8) mu = square(eta);              // sqrt link
    else if (link == 9) mu = inv(sqrt(eta));           // 1/mu^2 link
    else reject("Link not supported");
    // Compute likelihood
    if (dist == 1) {     // Bernoulli
      return dot_product(y, log(mu)) + dot_product(1 - y, log(1 - mu));
    }
    else if (dist == 2) {  // Poisson
      return dot_product(y, log(mu)) - sum(mu + lgamma(y + 1));
    }
    else if (dist == 3) {  // Normal
      return normal_lpdf(y | mu, sqrt(phi) );
    }
    else if (dist == 4) { // Gamma
      return gamma_lpdf(y | tau, tau * mu );
    }
    else if (dist == 5) { // Inverse-Gaussian
      return
          0.5 * (
              n * (log(tau) - 1.8378770664) - 3 * sum(log(y))
            - tau * dot_self( (y .* inv(mu) - 1) .* inv_sqrt(y) )
          );
    }
    else reject("Distribution not supported");
    return 0; // never reached;
  }
}
data {
  int<lower=0>         n;
  int<lower=0>         n0;
  int<lower=0>         p;
  vector[n]            y;
  matrix[n,p]          X;
  vector[n0]           y0;
  matrix[n0,p]         X0;
  vector[p]            beta0_mean;
  matrix[p,p]          beta0_cov;
  real<lower=0>        disp_shape;
  real<lower=0>        disp_scale;
  real<lower=0>        disp0_shape;
  real<lower=0>        disp0_scale;
  vector<lower=0>[p]   tau;
  int<lower=1,upper=5> dist;
  int<lower=1,upper=9> link;
  vector[n]            offset;
  vector[n0]           offset0;
}
transformed data{
  vector[p] comm_sd = inv_sqrt(tau);   // sd hyperparameter for commensurate prior
}
parameters {
  vector[p] beta;
  real<lower=0> dispersion[(dist > 2) ? 1 :  0];
  vector[p] beta0;
  real<lower=0> dispersion0[(dist > 2) ? 1 :  0];
}
model {
  beta0   ~ multi_normal(beta0_mean, beta0_cov);               // initial prior for beta
  beta    ~ normal(beta0, comm_sd);                            // commensurate prior for beta
  if ( dist <= 2 ) {
    target += glm_lp(y,  beta,  1.0, X,  dist, link, offset);    // current data likelihood
    target += glm_lp(y0, beta0, 1.0, X0, dist, link, offset0);   // historical data likelihood
  }
  else {
    dispersion0 ~ inv_gamma(disp0_shape, disp0_scale);                    // prior for historical dispersion
    dispersion  ~ inv_gamma(disp_shape,  disp_scale);                     // prior for dispersion
    target += glm_lp(y,  beta,  dispersion[1],  X,  dist, link, offset);  // current data likelihood
    target += glm_lp(y0, beta0, dispersion0[1], X0, dist, link, offset0); // historical data likelihood
  }
}

