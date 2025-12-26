data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  int y[N];
  vector[K] mu;
  matrix[K, K] Sigma; // (X'WX)^-1
  real<lower=0> g;    // Pass this as data, e.g., N
}
transformed data {
  matrix[K, K] L_prior;
  // Scale the covariance by g
  L_prior = cholesky_decompose(g * Sigma + 1e-8 * diag_matrix(rep_vector(1.0, K)));
}
parameters {
  vector[K] z;
}
transformed parameters {
  vector[K] beta = mu + L_prior * z;
}
model {
  z ~ std_normal();
  y ~ bernoulli_logit(X * beta);
}
