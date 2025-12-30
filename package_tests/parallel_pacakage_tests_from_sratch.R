pacman::p_load(clogitR, dplyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, doParallel, ggplot2, geepack, glmmTMB, rstan) #doParallel
options(error = recover)
rm(list = ls())
############### parameters ###############
num_cores = availableCores() - 6
Nsim = 100
external_nsim = 100000
ns = c(100, 250, 500)
beta_Ts = c(0,1)
X_styles = c("non-correlated") #"correlated", 
true_funtions = c("linear", "non-linear")
regress_on_Xs = c("one", "two") #"all", "one", "none"

#Bayesian_prior_for_betaTs = c(TRUE, FALSE)
sm = rstan::stan_model("mvn_logistic.stan")
sm_g = rstan::stan_model("mvn_logistic_gprior.stan")
sm_PMP = rstan::stan_model("mvn_logistic_PMP.stan")
sm_hybrid = rstan::stan_model("mvn_logistic_Hybrid.stan")

params = expand.grid(
  nsim = 1:Nsim,
  beta_T = beta_Ts,
  n = ns,
  X_style = X_styles
)

params = params %>%
  arrange(nsim, beta_T, n, X_style) 

############### CODE ###############

Bayesian_Clogit = function(y_dis, X_dis, w_dis, y_con, X_con, w_con, prior_type, concordant_fit) {
  
  if (concordant_fit == "reg") {
    fit_con = glm(y_con ~ w_con + X_con, family = "binomial")
    b_con = summary(fit_con)$coefficients[,1]
    Sigma_con = pmin(vcov(fit_con), 20)
    eps = 1e-6         #added for stabilibty
    Sigma_con = Sigma_con + diag(eps, nrow(Sigma_con))
    
    b_con = c(0, b_con[-c(1,2)])
    Sigma_con = Sigma_con[-1,-1]; Sigma_con[1,] = 0; Sigma_con[, 1] = 0; Sigma_con[1,1] = 20
  } else if (concordant_fit == "GEE") {
    strat_con = rep(1:(nrow(X_con)/2), each = 2)
    fit_con = geeglm(
      y_con ~ w_con + X_con,
      id    = strat_con,
      family = binomial(link = "logit"),
      corstr = "exchangeable",
      data   = data.frame(y_con, w_con, X_con, strat_con)
    )
    
    b_con = summary(fit_con)$coefficients[,1]
    Sigma_con = pmin(vcov(fit_con), 20)
    eps = 1e-6         #added for stabilibty
    Sigma_con = Sigma_con + diag(eps, nrow(Sigma_con))
    
    b_con = c(0, b_con[-c(1,2)])
    Sigma_con = Sigma_con[-1,-1]; Sigma_con[1,] = 0; Sigma_con[, 1] = 0; Sigma_con[1,1] = 20
  }
  
  y_dis_0_1 = ifelse(y_dis == -1, 0, 1)
  wX_dis = cbind(w_dis, X_dis)
  
  if (prior_type == "normal") {
    
    if (all(diag(Sigma_con) == 20)) {
      ret = list()
      ret$betaT = NA
      ret$ssq_beta_T = NA
      ret$reject = NA
      ret$pval = NA
      return(ret)#model blew up
    }
    
    data_list = list(
      N = nrow(wX_dis),
      K = ncol(wX_dis),
      X = wX_dis,
      y = y_dis_0_1,
      mu = b_con,        # example prior mean
      Sigma = Sigma_con    # example covariance (wide prior)
    )
    
    discordant_model = tryCatch({
      #for other aphas you would need to get the posterior for treatment effect and take the quarantines of that. you can get that from: post <- rstan::extract(fit); w_samples <- post$beta[, 1]
      summary(rstan::sampling(sm, data = data_list, refresh = 0, chains = 1))$summary["beta[1]", c("mean", "sd", "2.5%", "97.5%")]
    }, error = function(e) {
      warning(sprintf("stan_glm failed: %s", e$message))
      NULL
    })
    
  } else if (prior_type == "G prior") {
    
    if (all(diag(Sigma_con) == 20)) {
      ret = list()
      ret$betaT = NA
      ret$ssq_beta_T = NA
      ret$reject = NA
      ret$pval = NA
      return(ret)#model blew up
    }
    
    data_list = list(
      N = nrow(wX_dis),
      K = ncol(wX_dis),
      X = wX_dis,
      y = y_dis_0_1,
      mu = b_con,
      Sigma = Sigma_con
    )
    
    discordant_model = tryCatch({
      #for other aphas you would need to get the posterior for treatment effect and take the quarantines of that. you can get that from: post <- rstan::extract(fit); w_samples <- post$beta[, 1]
      summary(rstan::sampling(sm_g, data = data_list, refresh = 0, chains = 1))$summary["beta[1]", c("mean", "sd", "2.5%", "97.5%")]
    }, error = function(e) {
      warning(sprintf("stan_glm failed: %s", e$message))
      NULL
    })
    
  } else if (prior_type == "PMP") {
    
    # 3. Use the CONCORDANT covariance for the other covariates (X)
    # Exclude intercept [1] and treatment [2]
    mu_X_prior = matrix(b_con[-c(1)], ncol = 1)
    Sigma_X_prior = as.matrix(Sigma_con[-c(1), -c(1)])
    
    data_list = list(
      N = nrow(wX_dis),
      K = ncol(wX_dis),
      X = wX_dis,
      y = y_dis_0_1,
      mu_T = 0,  
      V_T = 20,
      mu_X = mu_X_prior, # Ensure vector for Stan
      Sigma_X = Sigma_X_prior
    )
    
    discordant_model = tryCatch({
      #for other aphas you would need to get the posterior for treatment effect and take the quarantines of that. you can get that from: post <- rstan::extract(fit); w_samples <- post$beta[, 1]
      summary(rstan::sampling(sm_PMP, data = data_list, refresh = 0, chains = 1))$summary["beta_T", c("mean", "sd", "2.5%", "97.5%")]
    }, error = function(e) {
      warning(sprintf("stan_glm failed: %s", e$message))
      NULL
    })
    
  } else if (prior_type == "Hybrid") {
    
    proj_matrix = X_dis %*% solve(t(X_dis) %*% X_dis) %*% t(X_dis)
    w_dis_ortho = w_dis - proj_matrix %*% w_dis
    
    
    # Prepare the data list for Stan
    data_list = list(
      N = nrow(X_dis),
      P = ncol(X_dis),
      y = y_dis_0_1,
      xw = as.vector(w_dis_ortho),
      X = X_dis,
      mu_A = as.array(b_con[-1]), 
      Sigma_A = as.matrix(Sigma_con[-1, -1, drop = FALSE])
    )
    
    # # 2. Build the data list to match the Stan model variables
    # data_list = list(
    #   N = nrow(wX_dis),
    #   K = ncol(wX_dis),
    #   X = wX_dis,
    #   y = y_dis_0_1,
    #   mu_T = 0,         # Center treatment prior at 0
    #   V_T = 20,        # High variance for treatment (non-informative)
    #   mu_X = mu_X_prior, 
    #   Sigma_X = Sigma_X_prior
    # )
    
    discordant_model = tryCatch({
      summary(rstan::sampling(sm_hybrid, data = data_list, refresh = 0, chains = 1))$summary["beta_w", c("mean", "sd", "2.5%", "97.5%")]
    }, error = function(e) {
      warning(sprintf("Hybrid PMP failed: %s", e$message))
      NULL
    })
  }
  
  ret = list(
    betaT = NA,
    ssq_beta_T = NA,
    reject = NA,
    pval = NA
  )
  
  if (!is.null(discordant_model)) {
    ret$betaT = discordant_model[1]
    ret$ssq_beta_T = discordant_model[2]
    ret$reject = !(discordant_model[3] < 0 & discordant_model[4] > 0)
    ret$pval = if (ret$reject) 0 else 1
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
  dis_idx =           matched_data$discordant_idx + 1
  
  if (regress_on_X %in% c("all", "one", "two")) {
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
    
    ########################### BAYESIAN ########################### 
    for (concordant_fit in c("reg", "GEE")) {
      beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA; pval_freq = NA
      if (discordant_viabele & concordant_viabele) {
        model = Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, "normal", concordant_fit)
        beta_hat_T = model$betaT; ssq_beta_hat_T = model$ssq_beta_T 
        pval = model$pval
        #pval_freq = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
      }
      res = rbind(res, data.frame(
        n = n,
        beta_T = beta_T,
        X_style = X_style,
        true_funtion = true_funtion,
        regress_on_X = regress_on_X,
        inference = paste0("bayesian_", concordant_fit),
        beta_hat_T = beta_hat_T,
        ssq_beta_hat_T = ssq_beta_hat_T,
        pval = pval
      ))
    }
    
    
    ########################### BAYESIAN G PRIOR ########################### 
    for (concordant_fit in c("reg", "GEE")) {
      beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA; pval_freq = NA
      if (discordant_viabele & concordant_viabele) {
        model = Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, "G prior", concordant_fit)
        beta_hat_T = model$betaT; ssq_beta_hat_T = model$ssq_beta_T 
        pval = model$pval
        #pval_freq = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
      }
      res = rbind(res, data.frame(
        n = n,
        beta_T = beta_T,
        X_style = X_style,
        true_funtion = true_funtion,
        regress_on_X = regress_on_X,
        inference = paste0("bayesian G prior_", concordant_fit),
        beta_hat_T = beta_hat_T,
        ssq_beta_hat_T = ssq_beta_hat_T,
        pval = pval
      ))
    }
    
    ########################### BAYESIAN PMP ########################### 
    # for (concordant_fit in c("reg", "GEE")) {
    #   beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA; pval_freq = NA
    #   if (discordant_viabele & concordant_viabele) {
    #     model = Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, "PMP", concordant_fit)
    #     beta_hat_T = model$betaT; ssq_beta_hat_T = model$ssq_beta_T 
    #     pval = model$pval
    #     #pval_freq = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    #   }
    #   res = rbind(res, data.frame(
    #     n = n,
    #     beta_T = beta_T,
    #     X_style = X_style,
    #     true_funtion = true_funtion,
    #     regress_on_X = regress_on_X,
    #     inference = paste0("bayesian PMP_", concordant_fit),
    #     beta_hat_T = beta_hat_T,
    #     ssq_beta_hat_T = ssq_beta_hat_T,
    #     pval = pval
    #   ))
    # }
    
    ########################### BAYESIAN Hybrid ########################### 
    for (concordant_fit in c("reg", "GEE")) {
      beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA; pval_freq = NA
      if (discordant_viabele & concordant_viabele) {
        model = Bayesian_Clogit(y_dis, X_dis, w_dis, y_con, X_con, w_con, "Hybrid", concordant_fit)
        beta_hat_T = model$betaT; ssq_beta_hat_T = model$ssq_beta_T 
        pval = model$pval
        #pval_freq = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
      }
      res = rbind(res, data.frame(
        n = n,
        beta_T = beta_T,
        X_style = X_style,
        true_funtion = true_funtion,
        regress_on_X = regress_on_X,
        inference = paste0("bayesian Hybrid_", concordant_fit),
        beta_hat_T = beta_hat_T,
        ssq_beta_hat_T = ssq_beta_hat_T,
        pval = pval
      ))
    }
    
    ########################### glmmTMB  ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    tryCatch({
      fit_tmb = glmmTMB(
        y ~ X + w + (1 | strat),
        family = binomial(),
        data   = data.frame(y, X, w, strat)
      )
      model = summary(fit_tmb)$coefficients$cond["w", c("Estimate", "Std. Error")]
      beta_hat_T = model[1]; ssq_beta_hat_T = model[2]
      pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }, error = function(e) {
      beta_hat_T <<- NA; ssq_beta_hat_T <<- NA; pval <<- NA
    })
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "glmmTMB",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
    
    ########################### GEE  ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    tryCatch({
      fit_gee = geeglm(
        y ~ X + w,
        id    = strat,
        family = binomial(link = "logit"),
        corstr = "exchangeable",
        data   = data.frame(y, X, w, strat)
      )
      model = summary(fit_gee)$coefficients["w", c("Estimate", "Std.err")]
      beta_hat_T = as.numeric(model[1]); ssq_beta_hat_T = as.numeric(model[2])
      pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }, error = function(e) {
      beta_hat_T <<- NA; ssq_beta_hat_T <<- NA; pval <<- NA
    })
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "GEE",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
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
    
    ########################### glmmTMB  ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    tryCatch({
      fit_tmb = glmmTMB(
        y ~ w + (1 | strat),
        family = binomial(),
        data   = data.frame(y, X, w, strat)
      )
      model = summary(fit_tmb)$coefficients$cond["w", c("Estimate", "Std. Error")]
      beta_hat_T = model[1]; ssq_beta_hat_T = model[2]
      pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }, error = function(e) {
      beta_hat_T <<- NA; ssq_beta_hat_T <<- NA; pval <<- NA
    })
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "glmmTMB",
      beta_hat_T = beta_hat_T,
      ssq_beta_hat_T = ssq_beta_hat_T,
      pval = pval
    ))
    
    ########################### GEE  ########################### 
    beta_hat_T = NA; ssq_beta_hat_T = NA; pval = NA
    tryCatch({
      fit_gee = geeglm(
        y ~ w,
        id    = strat,
        family = binomial(link = "logit"),
        corstr = "exchangeable",
        data   = data.frame(y, X, w, strat)
      )
      model = summary(fit_gee)$coefficients["w", c("Estimate", "Std.err")]
      beta_hat_T = as.numeric(model[1]); ssq_beta_hat_T = as.numeric(model[2])
      pval = 2 * pnorm(min(c(-1,1) * (beta_hat_T / ssq_beta_hat_T)))
    }, error = function(e) {
      beta_hat_T <<- NA; ssq_beta_hat_T <<- NA; pval <<- NA
    })
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      X_style = X_style,
      true_funtion = true_funtion,
      regress_on_X = regress_on_X,
      inference = "GEE",
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
  
  if (X_style == "correlated") { ############## correlated is the old style
    Sigma = 1 * (matrix(0.5, nrow = 6, ncol = 6) + diag(1 - 0.5, 6))
    X = MASS::mvrnorm(n/2, rep(0, 6), Sigma)
    X = pnorm(X)
    X = matrix(2*X - 1, ncol = 6)
  } else {
    X = matrix(runif((n/2) * 6, min = -1, max = 1), ncol = 6)
    X_plus_eps = X + matrix(rnorm((n/2) * 6, 0, 0.1),ncol = 6)
    combined = rbind(X, X_plus_eps)
    ids = order(c(1:(n/2), 1:(n/2)))
    X = combined[ids, ]
    X[,1] = runif(n, min = -1, max = 1)
    rm(X_plus_eps, combined, ids)
  }
  # df = data.frame(cbind(id = 1:n, X))
  # df.dist = gendistance(data.frame(df[, -1]), idcol = 1)
  # df.mdm = distancematrix(df.dist)
  # df.match = nonbimatch(df.mdm)
  # 
  # T_inx = df.match$halves[,2]
  # C_ind = df.match$halves[,4]
  # X = X[c(rbind(T_inx, C_ind)), ] #zip them together, so the mathces should be 1,1,2,2,3,3...
  w = c(rbind(replicate(n/2, sample(c(0, 1)), simplify = TRUE)))
  strat = rep(1:(n/2), each = 2)
  
  
  for (true_funtion in true_funtions) {
    if (true_funtion == "linear") {
      beta_X = c(1, 1, 1, 1, 1, 1)
      beta_0 = -0.5
      probs = 1 / (1 + exp(-(beta_0 + (as.matrix(X) %*% beta_X) + beta_T * w)))
    } else {
      f_x = sin(pi * X[, 1] * X[, 2]) + X[,3]^3 + X[, 4]^2 + X[, 5]^2
      probs = 1 / (1 + exp(-(f_x + beta_T * w)))
    }
    y = rbinom(n, 1, probs)
    
    # df = data.frame(y = y, probs = probs)
    # ggplot(df, aes(x = probs, fill = factor(y))) +
    #   geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
    #   labs(x = "Predicted probability", fill = "Outcome") +
    #   scale_fill_manual(values = c("0" = "red", "1" = "blue")) +
    #   theme_minimal()
    
    for (regress_on_X in regress_on_Xs) {
      if (regress_on_X == "one") {
        X_run = X[,1, drop = FALSE]
      } else if (regress_on_X == "two") {
        X_run = X[,c(1,2), drop = FALSE]
      } else {
        X_run = X
      }
      one_res = Do_Inference(y, X_run, w, strat, beta_T, n, X_style, true_funtion, regress_on_X)
      BIG_res = rbind(BIG_res, one_res)
    }
  }
  return(BIG_res)
}


for (j in 1:120) {
  cat("################", j, "################\n")
  beta_T = params[j,]$beta_T
  n = params[j,]$n
  X_style = params[j,]$X_style
  print(Run_sim(beta_T = beta_T, n = n, X_style = X_style)); cat('\n')
}

############### SIM SET UP ###############

handlers(global = TRUE)
handlers("txtprogressbar")

registerDoFuture()
plan(multisession, workers = num_cores)

############### SIM ###############

for (e_nsim in 22:external_nsim) {

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
  write.csv(results, file = paste0("C:/temp/clogitR_kap_test_from_scratch/", Nsim, "_", e_nsim, ".csv"), row.names = FALSE)
  rm(results); gc()
}

plan(sequential)


############### COMPILE RESULTS ###############

results = read.csv("C:/temp/clogitR_kap_test_from_scratch/100_1.csv")

sum = 1
for (i in 2:475) {
  file_path <- paste0("C:/temp/clogitR_kap_test_from_scratch/100_", i, ".csv")
  if (file.exists(file_path)) {
    sum = sum +1
    message("Reading file ", i)
    temp <- read.csv(file_path)
    results <- rbind(results, temp)
  } else {
    message("Skipping missing file ", i)
  }
}

# for (i in 2:475){
#   print(i)
#   results = rbind(results, read.csv(paste0("C:/temp/clogitR_kap_test_from_scratch/1000_", i, ".csv")))
# }

results$X = NULL

res_mod = results %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(beta_T, true_funtion, regress_on_X, n, inference, X_style) %>%
  summarize(
    num_na = sum(is.na(pval)), 
    num_real = sum(!is.na(pval)),
    mse = mean(sq_err, na.rm = TRUE),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    mean_beta_hat_T = mean(beta_hat_T, na.rm = TRUE),
    mean_sq_beta_hat_T = mean(sqrt(ssq_beta_hat_T), trim = 0.001, na.rm = TRUE),
    .groups = "drop")


write.csv(res_mod, file = "C:/temp/clogitR_kap_test_from_scratch/combined_4600.csv", row.names = FALSE)
