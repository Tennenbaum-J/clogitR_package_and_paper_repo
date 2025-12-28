data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  int<lower=0,upper=1> y[N];
  real mu_T;
  real<lower=0> V_T;
  matrix[K-1, 1] mu_X;      // Changed to matrix[K-1, 1] as requested
  matrix[K-1, K-1] Sigma_X;
}
parameters {
  real beta_T;
  vector[K-1] beta_X;
}
model {
  // 1. FREEDOM FUNCTION g(lambda)
  // Converting mu_X to a vector inside the call to match the parameter type
  if (K > 1) {
    beta_X ~ multi_normal(to_vector(mu_X), Sigma_X); //
  }
  beta_T ~ normal(mu_T, sqrt(V_T));

  // 2. THE PMP COMPONENT (Partial Fisher Information)
  vector[K] beta_full;
  beta_full[1] = beta_T;
  if (K > 1) {
    beta_full[2:K] = beta_X;
  }
  
  vector[N] p = inv_logit(X * beta_full);
  vector[N] w_vec = p .* (1.0 - p); // weights w_ii = p_i(1-p_i)
  
  matrix[K, K] Info;
  for (i in 1:K) {
    for (j in i:K) {
      Info[i, j] = dot_product(X[, i] .* w_vec, X[, j]); //
      Info[j, i] = Info[i, j]; 
    }
  }
  
  real partial_info;
  if (K > 1) {
    real I_psi_psi = Info[1, 1];
    vector[K-1] I_psi_lambda = Info[2:K, 1];
    matrix[K-1, K-1] I_lambda_lambda = Info[2:K, 2:K] + diag_matrix(rep_vector(1e-9, K-1));
    
    // PMP formula: I_psi_psi - I_psi_lambda * inv(I_lambda_lambda) * I_psi_lambda
    partial_info = I_psi_psi - dot_product(I_psi_lambda, mdivide_left_spd(I_lambda_lambda, I_psi_lambda));
  } else {
    partial_info = Info[1, 1];
  }
  
  if (partial_info > 1e-12) {
    target += 0.5 * log(partial_info); // The hybrid logic
  }

  // 3. LIKELIHOOD
  y ~ bernoulli_logit(X * beta_full);
}
