pacman::p_load(clogitR, dplyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, doParallel) #doParallel
#devtools::load_all("C:/Users/Jacob/clogitR_package_and_paper_repo/clogitR")
rm(list = ls())
set.seed(1986)
options(error = recover)

nsim_exact_test = 501
num_cores = availableCores()-10
external_nsim = 10000
Nsim = 500

beta_Ts = c(0, 1)
ns = c(100, 250, 500)

params = expand.grid(
  i = 1:Nsim,
  beta_T = beta_Ts,
  n = ns
)
params = params %>%
  arrange(i, beta_T, n) 

run_simulation = function(i, beta_T, n){ 
  #beta_T = 1; n = 100
  res = data.frame(
    n = numeric(),
    beta_T = numeric(),
    inference = character(),
    beta_hat_T = numeric(),
    ssq_hat_T = numeric(),
    pval = numeric()
  )

  probs = rep(NA, n)
  y = rep(NA, n)
  X = sort(rnorm(n, mean = 0.6, sd = 0.2))
  w = rep(c(0, 1), times = n / 2)
  strat = rep(1:(n/2), each = 2)
  beta_0 = -1
  beta_x = 1
  
  probs = 1 / (1 + exp(-(beta_0 + beta_x * X + beta_T * w)))
  y = rbinom(n, 1, probs)
  # cbind(y, w)
  
  matched_data =
    process_matched_pairs(
      strata = strat,
      y = y,
      X = data.matrix(X),
      treatment = w
    )

  reservoir_X =         matched_data$X_reservoir_concordant
  reservoir_y =         matched_data$y_reservoir_concordant
  reservoir_treatment = matched_data$treatment_reservoir_concordant
  diffs_X =             matched_data$X_diffs_discordant
  diffs_y =             matched_data$y_diffs_discordant
  diffs_treatment =     matched_data$treatment_diffs_discordant
  discordant_idx =      matched_data$discordant_idx
  
  # summary(glm(y[discordant_idx] ~ w[discordant_idx] + X[discordant_idx], family = "binomial"))
  # specil = 2 * y[discordant_idx] - 1
  # 
  # fast_conditional_logistic_regression_with_var_cpp(X_diff = cbind(1, w[discordant_idx], X[discordant_idx]), y_diff = specil, j = 2)
  
  discordant_no_diffs_model = 
    fastClogit(discordant_Xdiffs = NULL,
               discordant_ydiffs = NULL,
               discordnat_Treatmentdiffs = NULL,
               concordnat_X = cbind(1, X[discordant_idx]),
               concordnat_y = y[discordant_idx],
               concordnat_Treatment = w[discordant_idx])
  dnd_beta_hat = discordant_no_diffs_model$concordnat_betaT
  ssq_dnd_beta_hat = discordant_no_diffs_model$concordnat_ssq_b
  dnd_z_stat = c(-1,1) * (dnd_beta_hat / sqrt(ssq_dnd_beta_hat))
  dnd_prob = pnorm(dnd_z_stat)
  dnd_pval = 2 * min(dnd_prob)
  
  fastClogit(discordant_Xdiffs = diffs_X,
             discordant_ydiffs = diffs_y,
             discordnat_Treatmentdiffs = diffs_treatment,
             concordnat_X = NULL,
             concordnat_y = NULL,
             concordnat_Treatment = NULL)
  
  mixed_model =
    fastClogit(discordant_Xdiffs = diffs_X,
               discordant_ydiffs = diffs_y,
               discordnat_Treatmentdiffs = diffs_treatment,
               concordnat_X = reservoir_X,
               concordnat_y = reservoir_y,
               concordnat_Treatment = reservoir_treatment)
  m_beta_hat = mixed_model$mixed_betaT
  ssq_m_beta_hat = mixed_model$mixed_ssq_b
  m_z_stat = c(-1,1) * (m_beta_hat / sqrt(ssq_m_beta_hat))
  m_prob = pnorm(m_z_stat)
  m_pval = 2 * min(m_prob)
  
  d_beta_hat = mixed_model$discordnat_betaT
  ssq_d_beta_hat = mixed_model$discordnat_ssq_b
  d_z_stat = c(-1,1) * (d_beta_hat / sqrt(ssq_d_beta_hat))
  d_prob = pnorm(d_z_stat)
  d_pval = 2 * min(d_prob)
  
  c_beta_hat = mixed_model$concordnat_betaT
  ssq_c_beta_hat = mixed_model$concordnat_ssq_b
  c_z_stat = c(-1,1) * (c_beta_hat / sqrt(ssq_c_beta_hat))
  c_prob = pnorm(c_z_stat)
  c_pval = 2 * min(c_prob)
  
  logit_model = 
    fastClogit(discordant_Xdiffs = NULL,
               discordant_ydiffs = NULL,
               discordnat_Treatmentdiffs = NULL,
               concordnat_X = data.matrix(X),
               concordnat_y = y,
               concordnat_Treatment = w)
  l_beta_hat = logit_model$concordnat_betaT
  ssq_l_beta_hat = logit_model$concordnat_ssq_b
  l_z_stat = c(-1,1) * (l_beta_hat / sqrt(ssq_l_beta_hat))
  l_prob = pnorm(l_z_stat)
  l_pval = 2 * min(l_prob)
  
  
  bayesian_model =
    bayesianClogit(discordant_Xdiffs = diffs_X,
                   discordant_ydiffs = diffs_y,
                   discordnat_Treatmentdiffs = diffs_treatment,
                   concordnat_X = reservoir_X,
                   concordnat_y = reservoir_y,
                   concordnat_Treatment = reservoir_treatment)
  b_beta_hat = bayesian_model$discordnat_betaT
  ssq_b_beta_hat = bayesian_model$discordnat_ssq_b
  b_z_stat = c(-1,1) * (b_beta_hat / sqrt(ssq_b_beta_hat))
  b_prob = pnorm(b_z_stat)
  b_pval = 2 * min(b_prob)
  
  res = rbind(res, data.frame(
    n = n,
    beta_T = beta_T,
    infrence = "discordant no diffs",
    beta_hat_T = dnd_beta_hat,
    ssq_hat_T = ssq_dnd_beta_hat,
    pval = dnd_pval
  ))
  
  res = rbind(res, data.frame(
    n = n,
    beta_T = beta_T,
    infrence = "mixed",
    beta_hat_T = m_beta_hat,
    ssq_hat_T = ssq_m_beta_hat,
    pval = m_pval
  ))
  
  res = rbind(res, data.frame(
    n = n,
    beta_T = beta_T,
    infrence = "discordant",
    beta_hat_T = d_beta_hat,
    ssq_hat_T = ssq_d_beta_hat,
    pval = d_pval
  ))
  
  res = rbind(res, data.frame(
    n = n,
    beta_T = beta_T,
    infrence = "concordant",
    beta_hat_T = c_beta_hat,
    ssq_hat_T = ssq_c_beta_hat,
    pval = c_pval
  ))

  res = rbind(res, data.frame(
    n = n,
    beta_T = beta_T,
    infrence = "logit",
    beta_hat_T = l_beta_hat,
    ssq_hat_T = ssq_l_beta_hat,
    pval = l_pval
  ))


  res = rbind(res, data.frame(
    n = n,
    beta_T = beta_T,
    infrence = "bayesian",
    beta_hat_T = b_beta_hat,
    ssq_hat_T = ssq_b_beta_hat,
    pval = b_pval
  ))
  
  return(res)
}


handlers(global = TRUE)
handlers("txtprogressbar")


registerDoFuture()
plan(multisession, workers = num_cores)


start_time = Sys.time()

for (e_nsim in 1:external_nsim){
  with_progress({
  
    prog = progressor(along = 1:nrow(params))
  
    results = foreach(row = iter(params, by = "row"),
                      .combine = rbind,
                      .packages = c("clogitR", "nbpMatching", "data.table",
                                    "dplyr", "MASS", "Rcpp")) %dorng% {
      i = row$i
      n = row$n
      beta_T = row$beta_T
      
      #cat(glue::glue("Running i={i}"), '\n')
      res = tryCatch({
        out = run_simulation(i, beta_T, n)
        #cat("Successfully ran simulation")
        prog()
        out
      }, error = function(e) {
        #cat(glue::glue("Error in sim i={i}: {e$message}"), '\n')
        prog()  # still update progress bar even if it fails
        NULL    # return NULL if failed, will be dropped in rbind
      })
    }
  })
  write.csv(results, file = paste0("C:/temp/clogitR_kap_test_2/500_", e_nsim, ".csv"), row.names = FALSE)
  rm(results); gc()
}


plan(sequential)
end_time = Sys.time()


results = read.csv("C:/temp/clogitR_kap_test_2/500_1.csv")
for (i in 2:60){
  results = rbind(results, read.csv(paste0("C:/temp/clogitR_kap_test_2/500_", i, ".csv")))
}
results$X = NULL

res_mod = results %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(beta_T, infrence, n) %>%
  summarize(
    num_na = sum(is.na(pval)), 
    mse = mean(sq_err, na.rm = TRUE),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    mean_betaT_hat = mean(beta_hat_T, na.rm = TRUE),
    mean_sq_betaT_hat = mean(sqrt(ssq_hat_T), trim = 0.001, na.rm = TRUE),
    .groups = "drop")


write.csv(results, file = "C:/Users/Jacob/clogitR_package_and_paper_repo/package_tests/results.csv", row.names = FALSE)

summary(results[results$infrence == "discordant no diffs", "ssq_hat_T"])
res_mod[res_mod$infrence %in% c("discordant no diffs", "discordant"), c("beta_T", "infrence", "percent_reject", "mean_betaT_hat", "mean_sq_betaT_hat")]

results[(results$beta_T == 1) & 
          (results$betas == "0.5-2100") & 
          (results$rho == 0) & 
          (results$n == 100) & 
          (results$infrence == "matched"), ] |>
  (\(df) df[order(-df$beta_hat_T), ])()

library(ggplot2)

res_mod <- res_mod %>%
  mutate(beta_T = factor(beta_T))

# plot density distributions of mean_betaT_hat by inference type, faceted by beta_T
ggplot(results, aes(x = beta_T, color = infrence, fill = infrence)) +
  geom_histogram(alpha = 0.3) +
  facet_wrap(~ beta_T + n, scales = "free") +
  labs(
    title = "Density of Mean Beta_T Hat by Inference Type",
    x = expression(hat(beta)[T]),
    y = "Density",
    color = "Inference",
    fill = "Inference"
  ) +
  theme_minimal()

results_t = results[(results$infrence %in% c("discordant", "discordant no diffs")),]
results_c = results[(results$infrence %in% c("bayesian", "logit")),]

ggplot(results_t, aes(x = beta_hat_T, color = infrence, fill = infrence)) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.3, position = "identity", color = "black") +
  facet_wrap(~ beta_T + n, scales = "free") +
  labs(
    title = "Normalized Histogram of Beta_Hat_T by Inference Type",
    x = expression(hat(beta)[T]),
    y = "Density",
    color = "Inference",
    fill = "Inference"
  ) +
  xlim(c(-3,3))
  theme_minimal()

  ggplot(results_c, aes(x = beta_hat_T, color = infrence, fill = infrence)) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.3, position = "identity", color = "black") +
    facet_wrap(~ beta_T + n, scales = "free") +
    labs(
      title = "Normalized Histogram of Beta_Hat_T by Inference Type",
      x = expression(hat(beta)[T]),
      y = "Density",
      color = "Inference",
      fill = "Inference"
    ) +
    xlim(c(-3,3))
  theme_minimal()
