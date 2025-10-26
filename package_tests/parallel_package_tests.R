pacman::p_load(clogitR, dplyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, doParallel) #doParallel
#devtools::load_all("C:/Users/Jacob/clogitR_package_and_paper_repo/clogitR")
rm(list = ls())
set.seed(1986)
options(error = recover)

mu_x = 1
sigma_x = 1
sigma_e = 1
nsim_exact_test = 501
num_cores = availableCores()-1
external_nsim = 10000
Nsim = 50

beta_Ts = c(0, 1)
configs = list()
configs = c(configs, list(setting = list(p = 2, rho = 0,     betas = c(1, -3, 1.5, 0, 0))))
configs = c(configs, list(setting = list(p = 2, rho = 0.75,  betas = c(1, -3, 1.5, 0, 0))))
configs = c(configs, list(setting = list(p = 2, rho = 0,     betas = c(0.5, -2, 1, 0, 0))))
configs = c(configs, list(setting = list(p = 2, rho = 0.75,  betas = c(0.5, -2, 1, 0, 0))))
configs = c(configs, list(setting = list(p = 4, rho = 0,     betas = c(1, -2, 1, 1.5, -2))))
configs = c(configs, list(setting = list(p = 4, rho = 0.75,  betas = c(1, -2, 1, 1.5, -2))))
ns = c(250, 500, 1000)

params = expand.grid(
  i = 1:Nsim,
  beta_T = beta_Ts,
  config = configs,
  n = ns,
  KEEP.OUT.ATTRS = FALSE
)
params = params %>%
  arrange(i, beta_T, config, n)

run_simulation = function(i, beta_T, config, n){
  #i = 1; beta_T = 1; config = list(setting = list(p = 2, rho = 0, betas = c(1, -3, 1.5, 0, 0))); n = 500
  #cat("nsim:", i, ", 1\n")
  res = data.frame(
    i = numeric(),
    beta_T = numeric(),
    n = numeric(),
    betas = character(),
    rho = numeric(),
    inference = character(),
    beta_hat_T = numeric(),
    pval = numeric()
  )

  probs = rep(NA, n)
  y = rep(NA, n)
  errors = rnorm(n, 0, sigma_e)
  betas = config$setting$betas
  rho = config$setting$rho
  p = config$setting$p
  Sigma = sigma_x * (matrix(rho, nrow = p, ncol = p) + diag(1 - rho, p))
  X = data.table(MASS::mvrnorm(n, rep(mu_x, p), Sigma))

  df = data.frame(cbind(id = 1:n, X))
  df.dist = gendistance(df[, -1], idcol = 1)
  df.mdm = distancematrix(df.dist)
  df.match = nonbimatch(df.mdm)

  T_inx = df.match$halves[,2]
  C_inx = df.match$halves[,4]
  w_t = as.numeric(1:n %in% T_inx)

  strat = rep(0, n)

  for(pair_idx in seq_len(nrow(df.match$halves))) {
    strat[df.match$halves[pair_idx, 2]] = pair_idx
    strat[df.match$halves[pair_idx, 4]] = pair_idx
  }

  if (p == 2){
    z = betas[1] * X[, 1] +
      betas[2] * X[, 2] +
      betas[3] * X[, 1]^2 +
      betas[4] * X[, 2]^2 +
      betas[5] * X[, 1] * X[, 2]
  } else if (p == 4){
    z = betas[1] * X[, 1] +
      betas[2] * X[, 2] +
      betas[3] * X[, 1]^2 +
      betas[4] * X[, 3] +
      betas[5] * X[, 4]
  } else {
    z = X %*% betas
  }
  probs = 1 / (1 + exp(-(beta_T * w_t + as.numeric(unlist(z)) + errors)))
  y = rbinom(n, 1, probs)

  matched_data =
    process_matched_pairs(
      strata = strat,
      y = y,
      X = data.matrix(X),
      treatment = w_t
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
  ssq_c_beta_hat = mixed_model$concordnat_betaT
  c_z_stat = c(-1,1) * (c_beta_hat / sqrt(ssq_c_beta_hat))
  c_prob = pnorm(c_z_stat)
  c_pval = 2 * min(c_prob)
  
  logit_model = 
    fastClogit(discordant_Xdiffs = NULL,
               discordant_ydiffs = NULL,
               discordnat_Treatmentdiffs = NULL,
               concordnat_X = data.matrix(X),
               concordnat_y = y,
               concordnat_Treatment = w_t)
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
    betas = paste0(betas, collapse=""),
    beta_T = beta_T,
    n = n,
    rho = rho,
    infrence = "mixed",
    beta_hat_T = m_beta_hat,
    pval = m_pval
  ))
  
  res = rbind(res, data.frame(
    betas = paste0(betas, collapse=""),
    beta_T = beta_T,
    n = n,
    rho = rho,
    infrence = "discordant",
    beta_hat_T = d_beta_hat,
    pval = d_pval
  ))
  
  res = rbind(res, data.frame(
    betas = paste0(betas, collapse=""),
    beta_T = beta_T,
    n = n,
    rho = rho,
    infrence = "concordant",
    beta_hat_T = c_beta_hat,
    pval = c_pval
  ))

  res = rbind(res, data.frame(
    betas = paste0(betas, collapse=""),
    beta_T = beta_T,
    n = n,
    rho = rho,
    infrence = "logit",
    beta_hat_T = l_beta_hat,
    pval = l_pval
  ))


  res = rbind(res, data.frame(
    betas = paste0(betas, collapse=""),
    beta_T = beta_T,
    n = n,
    rho = rho,
    infrence = "bayesian",
    beta_hat_T = b_beta_hat,
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
      beta_T = row$beta_T
      config = row$config
      n = row$n
      cat(glue::glue("Running i={i}, n={n}"), '\n')
      res = tryCatch({
        out = run_simulation(i, beta_T, config, n)
        cat("Successfully ran simulation")
        prog()
        out
      }, error = function(e) {
        cat(glue::glue("Error in sim i={i}, n={n}: {e$message}"), '\n')
        prog()  # still update progress bar even if it fails
        NULL    # return NULL if failed, will be dropped in rbind
      })
    }
    write.csv(results, file = paste0("C:/temp/clogitR_power_2/2_50_", e_nsim, ".csv"), row.names = FALSE)
    rm(results); gc()
  }
})




plan(sequential)
end_time = Sys.time()


results = read.csv("C:/temp/clogitR_power_2/2_50_1.csv")
for (i in 2:585){
  results = rbind(results, read.csv(paste0("C:/temp/clogitR_power_2/2_50_", i, ".csv")))
}
results$X = NULL

res_mod = results %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(betas, rho, beta_T, n, infrence) %>%
  summarize(
    num_na = sum(is.na(pval)), 
    mse = mean(sq_err, na.rm = TRUE),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    .groups = "drop")


write.csv(res_mod, file = "C:/temp/clogitR_power_2/2_combined_29250.csv", row.names = FALSE)



results[(results$beta_T == 1) & 
          (results$betas == "0.5-2100") & 
          (results$rho == 0) & 
          (results$n == 100) & 
          (results$infrence == "matched"), ] |>
  (\(df) df[order(-df$beta_hat_T), ])()
