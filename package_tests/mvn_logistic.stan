data {
  int<lower=1> N;                      // number of observations
  int<lower=1> K;                      // number of predictors
  matrix[N, K] X;                      // design matrix
  int<lower=0,upper=1> y[N];           // binary outcomes

  vector[K] mu;                        // prior mean vector
  matrix[K, K] Sigma;                  // prior covariance matrix
}

transformed data {
  matrix[K, K] Sigma_chol;             // Cholesky factor of covariance
  Sigma_chol = cholesky_decompose(Sigma);
}

parameters {
  vector[K] beta;                      // regression coefficients
}

model {
  // Multivariate normal prior
  beta ~ multi_normal_cholesky(mu, Sigma_chol);

  // Logistic regression likelihood
  y ~ bernoulli_logit(X * beta);
}
