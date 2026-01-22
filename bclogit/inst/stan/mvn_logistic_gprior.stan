data {
  int<lower=1> N;                      // number of observations
  int<lower=1> K;                      // number of predictors
  matrix[N, K] X;                      // design matrix
  int<lower=0,upper=1> y[N];           // binary outcomes
  vector[K] mu;                        // prior mean vector
  matrix[K, K] Sigma;                  // prior covariance matrix
}

transformed data {
  matrix[K, K] L_Sigma = cholesky_decompose(Sigma);
}

parameters {
vector[K] z;                           // Helper: Standard Normal(0,1)
  real<lower=0> g;                     // The hyper-parameter g (the scaling factor)
}

transformed parameters {
  // Manual re-centering: beta = mu + sqrt(g) * L * z
  // This is mathematically the same as beta ~ N(mu, g*Sigma)
  // The reason we do this is to decouple g from beta
  vector[K] beta = mu + sqrt(g) * (L_Sigma * z);
}

model {
  // Hyper-prior: Zellner-Siow
  g ~ inv_gamma(0.5, N/2.0);
  // Prior for the helper (leads to the g-prior for beta)
  z ~ std_normal();
  // Likelihood
  y ~ bernoulli_logit(X * beta);
}
