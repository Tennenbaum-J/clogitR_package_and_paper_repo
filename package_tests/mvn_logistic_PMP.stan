data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  int<lower=0,upper=1> y[N];

  real mu_T;
  real<lower=0> V_T;

  matrix[K-1, 1] mu_X;
  matrix[K-1, K-1] Sigma_X;
}

parameters {
  real beta_T;
  vector[K-1] beta_X;
}

model {
  beta_T ~ normal(mu_T, sqrt(V_T));
  beta_X ~ multi_normal(to_vector(mu_X), Sigma_X);

  y ~ bernoulli_logit(
        beta_T * X[,1]
        + X[,2:K] * beta_X
      );
}
