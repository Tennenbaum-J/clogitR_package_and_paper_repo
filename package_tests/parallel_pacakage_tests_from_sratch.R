pacman::p_load(clogitR, dplyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, doParallel, ggplot2) #doParallel
options(error = recover)
rm(list = ls())
############### parameters ###############
num_cores = availableCores() - 4
Nsim = 1000
external_nsim = 100000
ns = c(100, 250, 500)
beta_Ts = c(0,1)
X_styles = c("correlated", "non-correlated")
true_funtions = c("linear", "non-linear")
regress_on_Xs = c("all", "one", "none")

#Bayesian_prior_for_betaTs = c(TRUE, FALSE)


params = expand.grid(
  nsim = 1:Nsim,
  beta_T = beta_Ts,
  n = ns,
  X_style = X_styles
)

params = params %>%
  arrange(nsim, beta_T, n, X_style) 

############### CODE ###############

Bayesian_Clogit = function(y_dis, X_dis, w_dis, y_con, X_con, w_con, inculde_prior_for_beta_T) {
  concordant_model = summary(glm(y_con ~ w_con + X_con, family = "binomial"))$coefficients[,c(1,2)]
  
  
  if (inculde_prior_for_beta_T) {
    prior_means = concordant_model[-1,1]
    prior_cov_no_intercept = pmin(concordant_model[-1,2], 20)
  } else {
    prior_means = c(0, concordant_model[-c(1,2),1])
    prior_cov_no_intercept = pmin(c(20, concordant_model[-c(1,2),2]), 20)
  }
  
  if (all(prior_cov_no_intercept == 20)) {
    ret = list()
    ret$betaT = NA
    ret$ssq_beta_T = NA
    ret$pval = NA
    return(ret)#model blew up
  }
  
  y_dis_0_1 = ifelse(y_dis == -1, 0, 1)
  discordant_model = tryCatch({
    rstanarm::stan_glm(
      y_dis_0_1 ~ 0 + w_dis + X_dis,
      family = binomial(link = "logit"),
      prior = rstanarm::normal(location = prior_means, scale = prior_cov_no_intercept),
      data = data.frame(y_dis_0_1, w_dis, X_dis),
      refresh = 0
    )
  }, error = function(e) {
    warning(sprintf("stan_glm failed: %s", e$message))
    NULL
  })
  
  ret = list()
  if (!is.null(discordant_model)){
    post = as.matrix(discordant_model)
    coef_name = grep("^w_dis", colnames(post), value = TRUE)
    w_samples = post[, coef_name]
    
    ret$betaT = mean(w_samples)
    ret$ssq_beta_T = sd(w_samples)
    ret$pval = 2 * min(mean(w_samples > 0), mean(w_samples < 0))
  }
  
  return(ret)
}

Do_Inference = function(y, X, w, strat, beta_T, n, X_style, true_funtion, regress_on_X) {
  res = data.frame(
    n = numeric(),
    beta_T = numeric(),
    X_style = character(),
    true_funtion = character(),
    regress_on_X = character(),
    inference = character(),
    beta_hat_T = numeric(),
    pval = numeric()
  )
  matched_data =
    process_matched_pairs(
      strata = strat,
      y = y,
      X = data.matrix(X),
      treatment = w
    )
  
  X_con =         matched_data$X_reservoir_concordant
  y_con =         matched_data$y_reservoir_concordant
  w_con = matched_data$treatment_reservoir_concordant
  X_dis =             matched_data$X_diffs_discordant
  y_dis =             matched_data$y_diffs_discordant
  w_dis =     matched_data$treatment_diffs_discordant
  dis_idx =      matched_data$discordant_idx
  
  if (regress_on_X %in% c("all", "one")) {
    discordant_viabele = if(length(y_dis) > ncol(X) + 7) { TRUE } else { FALSE }
    concordant_viabele = if(length(y_con) > ncol(X) + 7) { TRUE } else { FALSE }
    
    ########################### CLOGIT  ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    if (discordant_viabele) {
      y_dis_0_1 = ifelse(y_dis == -1, 0, 1)
      model = summary(glm(y_dis_0_1 ~ 0 + w_dis + X_dis, family = "binomial"))$coefficients[1,c(1,2)]
      beta_hat_T = model[1]; ssq_beta_hat_T = model[2]
      pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "clogit",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
    
    # ########################### DISCORDANT  ########################### 
    # beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    # if (discordant_viabele) {
    #   model = summary(glm(y[dis_idx] ~ w[dis_idx] + X[dis_idx,], family = "binomial"))$coefficients[2,c(1,2)]
    #   beta_hat_T = model[1]; ssq_beta_hat_T = model[2]
    #   pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    # }
    # res = rbind(res, data.frame(
    #   n = n,
    #   beta_T = beta_T,
    #   X_style = X_style,
    #   true_funtion = true_funtion,
    #   regress_on_X = regress_on_X,
    #   inference = "discordant",
    #   beta_hat_T = beta_hat_T,
    #   ssq_beta_hat_T = ssq_beta_hat_T,
    #   pval = pval
    # ))
    
    ########################### LOGIT  ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    if (TRUE) {
      model = summary(glm(y ~ w + X, family = "binomial"))$coefficients[2,c(1,2)]
      beta_hat_T = model[1]; ssq_beta_hat_T = model[2]
      pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "logit",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
    
    ########################### BAYESIAN NO T PRIOR ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA; pval_freq = NA
    if (discordant_viabele & concordant_viabele) {
      model = Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, FALSE)
      beta_hat_T = model$betaT; ssq_beta_hat_T = model$ssq_beta_T 
      pval = model$pval
      pval_freq = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "bayesian no T prior",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "bayesian no T prior pavl-freq",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval_freq
    ))
    
    ########################### BAYESIAN ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA; pval_freq = NA
    if (discordant_viabele & concordant_viabele) {
      model = Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, TRUE)
      beta_hat_T = model$betaT; ssq_beta_hat_T = model$ssq_beta_T 
      pval = model$pval
      pval_freq = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "bayesian",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
    
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "bayesian pavl-freq",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval_freq
    ))
    
  } else { #if no x then remove the x parameter
    discordant_viabele = if(length(y_dis) > 5) { TRUE } else { FALSE }
    
    ########################### CLOGIT  ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    if (discordant_viabele) {
      y_dis_0_1 = ifelse(y_dis == -1, 0, 1)
      model = summary(glm(y_dis_0_1 ~ 0 + w_dis, family = "binomial"))$coefficients[1,c(1,2)]
      beta_hat_T = model[1]; ssq_beta_hat_T = model[2]
      pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "clogit",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
    
    # ########################### DISCORDANT  ########################### 
    # beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    # if (discordant_viabele) {
    #   model = summary(glm(y[dis_idx] ~ w[dis_idx], family = "binomial"))$coefficients[2,c(1,2)]
    #   beta_hat_T = model[1]; ssq_beta_hat_T = model[2]
    #   pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    # }
    # res = rbind(res, data.frame(
    #   n = n,
    #   beta_T = beta_T,
    #   X_style = X_style,
    #   true_funtion = true_funtion,
    #   regress_on_X = regress_on_X,
    #   inference = "discordant",
    #   beta_hat_T = beta_hat_T,
    #   ssq_beta_hat_T = ssq_beta_hat_T,
    #   pval = pval
    # ))
    
    ########################### LOGIT  ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    if (discordant_viabele) {
      model = summary(glm(y ~ w, family = "binomial"))$coefficients[2,c(1,2)]
      beta_hat_T = model[1]; ssq_beta_hat_T = model[2]
      pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "logit",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
    
  }
  
  rownames(res) = NULL
  return(res)
}

Run_sim = function(beta_T, n, X_style) {
  BIG_res = data.frame()
  y = array(NA, n)
  probs = array(NA, n)
  
  if (X_style == "correlated") {
    Sigma = 1 * (matrix(0.5, nrow = 6, ncol = 6) + diag(1 - 0.5, 6))
    X = MASS::mvrnorm(n, rep(0, 6), Sigma) #error found in this line, the mean was set to 1, but it should have been 0
    X = pnorm(X)
    X = matrix(2*X - 1, ncol = 6)
  } else {
    X = matrix(runif(n * 6, min = -1, max = 1), ncol = 6)
  }
  df = data.frame(cbind(id = 1:n, X))
  df.dist = gendistance(data.frame(df[, -1]), idcol = 1)
  df.mdm = distancematrix(df.dist)
  df.match = nonbimatch(df.mdm)
  
  T_inx = df.match$halves[,2]
  C_ind = df.match$halves[,4]
  X = X[c(rbind(T_inx, C_ind)), ] #zip them together, so the mathces should be 1,1,2,2,3,3...
  w = c(rbind(replicate(n/2, sample(c(0, 1)), simplify = TRUE)))
  strat = rep(1:(n/2), each = 2)
  
  
  for (true_funtion in true_funtions) {
    if (true_funtion == "linear") {
      beta_X = c(3, 3, 3, 3, 3, 3)
      beta_0 = -1
      probs = 1 / (1 + exp(-(beta_0 + (as.matrix(X) %*% beta_X) + beta_T * w)))
    } else {
      f_x = sin(pi * X[, 1] * X[, 2]) + X[,3]^3 + X[, 4]^2 + X[, 5]^2
      probs = 1 / (1 + exp(-(f_x + beta_T * w)))
    }
    y = rbinom(n, 1, probs)
    
    # ggplot(df, aes(x = probs, fill = factor(y))) +
    #   geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
    #   labs(x = "Predicted probability", fill = "Outcome") +
    #   scale_fill_manual(values = c("0" = "red", "1" = "blue")) +
    #   theme_minimal()
    
    for (regress_on_X in regress_on_Xs) {
      if (regress_on_X == "one") {
        X_run = X[,1, drop = FALSE]
      } else {
        X_run = X
      }
      one_res = Do_Inference(y, X_run, w, strat, beta_T, n, X_style, true_funtion, regress_on_X)
      BIG_res = rbind(BIG_res, one_res)
    }
  }
  return(BIG_res)
}


# for (j in 1:1000) {
#   cat("################", j, "################\n")
#   beta_T = params[j,]$beta_T
#   n = params[j,]$n
#   X_style = params[j,]$X_style
#   print(Run_sim(beta_T = beta_T, n = n, X_style = X_style)); cat('\n')
# }

############### SIM SET UP ###############

handlers(global = TRUE)
handlers("txtprogressbar")

registerDoFuture()
plan(multisession, workers = num_cores)

############### SIM ###############

for (e_nsim in 1:external_nsim) {

  with_progress({
    prog = progressor(along = 1:nrow(params))

    results = foreach(row = iter(params, by = "row"),
                      .combine = rbind,
                      .packages = c("clogitR", "nbpMatching", "data.table",
                                    "dplyr", "MASS", "Rcpp", "rstanarm")) %dorng% {

      # extract parameters
      nsim = row$nsim
      beta_T = row$beta_T
      n = row$n
      X_style = row$X_style
      res = tryCatch({
        out = Run_sim(beta_T, n, X_style)
        #cat("Successfully ran simulation")
        prog()
        out
      }, error = function(e) {
        cat(glue::glue("Error in nsim={nsim}: {e$message}"), '\n')
        prog()  # still update progress bar even if it fails
        NULL    # return NULL if failed, will be dropped in rbind
      })
    }
  })
  write.csv(results, file = paste0("C:/temp/clogitR_kap_test_from_scratch/1000_", e_nsim, ".csv"), row.names = FALSE)
  rm(results); gc()
}

plan(sequential)


############### COMPILE RESULTS ###############

results = read.csv("C:/temp/clogitR_kap_test_from_scratch/1000_1.csv")

sum = 1
for (i in 2:475) {
  file_path <- paste0("C:/temp/clogitR_kap_test_from_scratch/1000_", i, ".csv")
  if (file.exists(file_path)) {
    sum = sum +1
    message("Reading file ", i)
    temp <- read.csv(file_path)
    results <- rbind(results, temp)
  } else {
    message("Skipping missing file ", i)
  }
}

for (i in 2:475){
  print(i)
  results = rbind(results, read.csv(paste0("C:/temp/clogitR_kap_test_from_scratch/1000_", i, ".csv")))
}

results$X = NULL

res_mod = results %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(beta_T, true_funtion, regress_on_X, n, inference) %>%
  summarize(
    num_na = sum(is.na(pval)), 
    num_real = sum(!is.na(pval)),
    mse = mean(sq_err, na.rm = TRUE),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    mean_beta_hat_T = mean(beta_hat_T, na.rm = TRUE),
    mean_sq_beta_hat_T = mean(sqrt(ssq_beta_hat_T), trim = 0.001, na.rm = TRUE),
    .groups = "drop")


write.csv(res_mod, file = "C:/temp/clogitR_kap_test_from_scratch/combined_13000.csv", row.names = FALSE)

