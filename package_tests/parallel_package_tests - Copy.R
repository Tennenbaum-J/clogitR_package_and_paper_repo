pacman::p_load(clogitR, dplyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, doParallel) #doParallel
#devtools::load_all("C:/Users/Jacob/clogitR_package_and_paper_repo/clogitR")
rm(list = ls())
set.seed(1986)
options(error = recover)

nsim_exact_test = 501
num_cores = availableCores()-1
external_nsim = 10000
Nsim = 500

beta_Ts = c(0, 1)
ns = c(100)

params = expand.grid(
  i = 1:Nsim,
  beta_T = beta_Ts,
  n = ns
)
params = params %>%
  arrange(i, beta_T, n) 

run_simulation = function(i, beta_T, n){ 
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
with_progress({

  for (e_nsim in 1:external_nsim){
    prog = progressor(along = 1:(nrow(params) * external_nsim))

    results = foreach(row = iter(params, by = "row"), .combine = rbind, .packages = c("clogitR", "nbpMatching", "data.table", "dplyr", "MASS", "Rcpp")) %dorng% {
      i = row$i
      n = row$n
      beta_T = row$beta_T
      cat(glue::glue("Running i={i}"), '\n')
      res = tryCatch({
        out = run_simulation(i, beta_T, n)
        cat("Successfully ran simulation")
        prog()
        out
      }, error = function(e) {
        cat(glue::glue("Error in sim i={i}: {e$message}"), '\n')
        prog()  # still update progress bar even if it fails
        NULL    # return NULL if failed, will be dropped in rbind
      })
    }
    write.csv(results, file = paste0("C:/temp/clogitR_kap_test/500_", e_nsim, ".csv"), row.names = FALSE)
    rm(results); gc()
  }
})




plan(sequential)
end_time = Sys.time()


results = read.csv("C:/temp/clogitR_kap_test/500_1.csv")
for (i in 2:120){
  results = rbind(results, read.csv(paste0("C:/temp/clogitR_kap_test/500_", i, ".csv")))
}
results$X = NULL

res_mod = results %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(beta_T, infrence) %>%
  summarize(
    num_na = sum(is.na(pval)), 
    mse = mean(sq_err, na.rm = TRUE),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    mean_betaT_hat = mean(beta_hat_T, na.rm = TRUE),
    mean_ssq_betaT_hat = mean(ssq_hat_T, trim = 0.0001, na.rm = TRUE),
    .groups = "drop")


write.csv(res_mod, file = "C:/temp/clogitR_kap_test/combined_60000.csv", row.names = FALSE)



results[(results$beta_T == 1) & 
          (results$betas == "0.5-2100") & 
          (results$rho == 0) & 
          (results$n == 100) & 
          (results$infrence == "matched"), ] |>
  (\(df) df[order(-df$beta_hat_T), ])()





