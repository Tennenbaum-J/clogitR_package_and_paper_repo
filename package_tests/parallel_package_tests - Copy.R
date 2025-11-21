pacman::p_load(clogitR, dplyr, data.table, doFuture, future, doRNG, foreach, progressr, doParallel, nbpMatching, doParallel, ggplot2) #doParallel
#devtools::load_all("C:/Users/Jacob/clogitR_package_and_paper_repo/clogitR")
rm(list = ls())
set.seed(1986)
options(error = recover)
#options(mc.cores = 1)

nsim_exact_test = 501
num_cores = availableCores()-5
external_nsim = 10000
Nsim = 100

beta_Ts = c(0, 1)
ns = c(100, 250, 500, 1000)
match_on_mores = c("only one", "corelated linear", "corelated non-linear", "not coralated linear", "not coralated non-linear")
regress_on_Xs = c("all", "one", "none")

params = expand.grid(
  i = 1:Nsim,
  beta_T = beta_Ts,
  n = ns,
  match_on_more = match_on_mores,
  regress_on_X = regress_on_Xs
)
params = params %>%
  arrange(i, beta_T, n, match_on_more, regress_on_X) 

params = params[-which(params$match_on_more == "only one" & params$regress_on_X == "all"),]

run_simulation = function(i, beta_T, n, match_on_more, regress_on_X){ 
  #beta_T = 1; n = 1000; match_on_more= "coralated"; regress_on_X = "all"
  probs = rep(NA, n)
  y = rep(NA, n)
  
  if (match_on_more == "only one") {
    X = matrix(sort(runif(n, -1,1)))
    w = rep(c(0, 1), times = n / 2)
    strat = rep(1:(n/2), each = 2)
    beta_0 = -1
    beta_x = 1
    
    probs = 1 / (1 + exp(-(beta_0 + beta_x * X + beta_T * w)))
    y = rbinom(n, 1, probs)
  } else {
    if (match_on_more %in% c("corelated linear", "corelated non-linear")) {
      Sigma = 1 * (matrix(0.5, nrow = 6, ncol = 6) + diag(1 - 0.5, 6))
      X = MASS::mvrnorm(n, rep(1, 6), Sigma)
      X = pnorm(X)
      X = matrix(2*X - 1, ncol = 6)
    } else {
      X = matrix(runif(n * 6, -1, 1), ncol = 6)
    }
    df = data.frame(cbind(id = 1:n, X))
    df.dist = gendistance(data.frame(df[, -1]), idcol = 1)
    df.mdm = distancematrix(df.dist)
    df.match = nonbimatch(df.mdm)
    
    T_inx = df.match$halves[,2]
    C_ind = df.match$halves[,4]
    X = X[c(rbind(T_inx, C_ind)), ]
    w = c(rbind(replicate(n/2, sample(c(0, 1)), simplify = TRUE)))
    strat = rep(1:(n/2), each = 2)
    
    
    if (match_on_more %in% c("corelated linear", "not coralated linear")) {
      beta_0 = -1
      beta_x = c(3, 3, 3, 3, 3, 3)
      probs = 1 / (1 + exp(-(beta_0 + (as.matrix(X) %*% beta_x) + beta_T * w)))
    } else {
      f_x = sin(pi * X[, 1] * X[, 2]) + X[,3]^3 + X[, 4]^2 + X[, 5]^2
      probs = 1 / (1 + exp(-(f_x + beta_T * w)))
    }
    y = rbinom(n, 1, probs)
    
    
    if (regress_on_X == "one") {
      X = X[,1]
    }
    if (FALSE) {
      #Nsim = 10000
      #pval_x = matrix(NA, Nsim, 1)
      #pval_w_x = array(NA, Nsim)
      #val_w = array(NA, Nsim)
      
      # for (nsim in 1 : Nsim) {
      #   w = c(rbind(replicate(n/2, sample(c(0, 1)), simplify = TRUE)))
      #   beta_0 = -1
      #   beta_x = c(3, 3, 3, 3, 3, 3)
      #   probs = 1 / (1 + exp(-(beta_0 + (as.matrix(X) %*% beta_x) + beta_T * w)))
      #   y = rbinom(n, 1, probs)
      #   
      #   obj = glm(y ~ data.matrix(X[,1]), family = "binomial")
      #   pval_x[nsim,] = coef(summary(obj))[2,4]
      #   obj = glm(y ~ w + data.matrix(X[,1]), family = "binomial")
      #   pval_w_x[nsim] = coef(summary(obj))[2,4]
      #   obj = glm(y ~ w, family = "binomial")
      #   pval_w[nsim] = coef(summary(obj))[2,4]
      # }
      # 
      # colMeans(pval_x < alpha)
      # mean(pval_w_x[!is.nan(pval_w_x)] < alpha)
      # mean(pval_w < alpha)
      # 
      # w = c(rbind(replicate(n/2, sample(c(0, 1)), simplify = TRUE)))
      # beta_0 = -1
      # beta_x = c(1, 2, -2, -3, 3, 0)
      # probs = 1 / (1 + exp(-(beta_0 + (as.matrix(X) %*% beta_x) + beta_T * w)))
      # y = rbinom(n, 1, probs)
      # 
      # df = data.frame(
      #   p = probs,
      #   y = y
      # )
      # 
      # ggplot(df, aes(x = p, fill = factor(y))) +
      #   geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
      #   labs(x = "Predicted probability", fill = "Outcome") +
      #   scale_fill_manual(values = c("0" = "red", "1" = "blue")) +
      #   theme_minimal()
      
      #drop the X to one col:
      #X = X[,1]
    }
  }
  
  return(computeAllData(strat, y, X, w, match_on_more, regress_on_X, n, beta_T))
}

computeAllData = function(strat, y, X, w, match_on_more, regress_on_X, n, beta_T) {
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
  if(regress_on_X == "all") {
    p = ncol(X)
  } else {
    p = 2
  }
  
  if (length(reservoir_y) < 2*p+5 || length(diffs_y) < 2*p+5){
    return(NULL)
  }
  
  res = data.frame(
    n = numeric(),
    beta_T = numeric(),
    match_on_more = character(),
    regress_on_X = character(),
    inference = character(),
    beta_hat_T = numeric(),
    ssq_hat_T = numeric(),
    pval = numeric()
  )
  
  if (regress_on_X %in% c("all", "one")) {
    discordant_no_diffs_model = 
      fastClogit(discordant_Xdiffs = NULL,
                 discordant_ydiffs = NULL,
                 discordnat_Treatmentdiffs = NULL,
                 concordnat_X = cbind(1, X[discordant_idx]),
                 concordnat_y = y[discordant_idx],
                 concordnat_Treatment = w[discordant_idx])
    dnd_beta_hat = discordant_no_diffs_model$concordnat_betaT
    ssq_dnd_beta_hat = discordant_no_diffs_model$concordnat_ssq_b
    dnd_pval = 2 * pnorm(min(c(-1,1) * (dnd_beta_hat / sqrt(ssq_dnd_beta_hat))))
    
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      match_on_more = match_on_more,
      regress_on_X = regress_on_X,
      infrence = "discordant no diffs",
      beta_hat_T = dnd_beta_hat,
      ssq_hat_T = ssq_dnd_beta_hat,
      pval = dnd_pval
    ))
    
    discordant_model =
      fastClogit(discordant_Xdiffs = diffs_X,
                 discordant_ydiffs = diffs_y,
                 discordnat_Treatmentdiffs = diffs_treatment,
                 concordnat_X = NULL,
                 concordnat_y = NULL,
                 concordnat_Treatment = NULL)
    
    
    d_beta_hat = discordant_model$discordnat_betaT
    ssq_d_beta_hat = discordant_model$discordnat_ssq_b
    d_pval = 2 * pnorm(min(c(-1,1) * (d_beta_hat / sqrt(ssq_d_beta_hat))))
    
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      match_on_more = match_on_more,
      regress_on_X = regress_on_X,
      infrence = "conditional logit",
      beta_hat_T = d_beta_hat,
      ssq_hat_T = ssq_d_beta_hat,
      pval = d_pval
    ))
    
    logit_model = 
      fastClogit(discordant_Xdiffs = NULL,
                 discordant_ydiffs = NULL,
                 discordnat_Treatmentdiffs = NULL,
                 concordnat_X = data.matrix(X),
                 concordnat_y = y,
                 concordnat_Treatment = w)
    l_beta_hat = logit_model$concordnat_betaT
    ssq_l_beta_hat = logit_model$concordnat_ssq_b
    l_pval = 2 * pnorm(min(c(-1,1) * (l_beta_hat / sqrt(ssq_l_beta_hat))))
    
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      match_on_more = match_on_more,
      regress_on_X = regress_on_X,
      infrence = "logit",
      beta_hat_T = l_beta_hat,
      ssq_hat_T = ssq_l_beta_hat,
      pval = l_pval
    ))
    
    bayesian_model =
      bayesianClogit(discordant_Xdiffs = diffs_X,
                     discordant_ydiffs = diffs_y,
                     discordnat_Treatmentdiffs = diffs_treatment,
                     concordnat_X = reservoir_X,
                     concordnat_y = reservoir_y,
                     concordnat_Treatment = reservoir_treatment)
    # Check if bayesian model actually produced results
    if (!is.null(bayesian_model$discordnat_betaT) && !is.null(bayesian_model$discordnat_ssq_b)) {
      b_beta_hat = bayesian_model$discordnat_betaT
      ssq_b_beta_hat = bayesian_model$discordnat_ssq_b
      b_pval = 2 * pnorm(min(c(-1,1) * (b_beta_hat / sqrt(ssq_b_beta_hat))))
      
      res = rbind(res, data.frame(
        n = n,
        beta_T = beta_T,
        match_on_more = match_on_more,
        regress_on_X = regress_on_X,
        infrence = "bayesian",
        beta_hat_T = b_beta_hat,
        ssq_hat_T = ssq_b_beta_hat,
        pval = b_pval
      ))
    } else {
      # Bayesian model failed, add NA row or skip
      res = rbind(res, data.frame(
        n = n,
        beta_T = beta_T,
        match_on_more = match_on_more,
        regress_on_X = regress_on_X,
        infrence = "bayesian",
        beta_hat_T = NA,
        ssq_hat_T = NA,
        pval = NA
      ))
    }
  } else {
    discordant_no_diffs_model_glm_no_X = summary(glm(y[discordant_idx] ~ w[discordant_idx], family = "binomial"))$coefficients[2,c(1,2)]
    dndgnX_beta_hat = discordant_no_diffs_model_glm_no_X[1]
    ssq_dndgnX_beta_hat = discordant_no_diffs_model_glm_no_X[2]^2
    dndgnX_pval = 2 * pnorm(min(c(-1,1) * (dndgnX_beta_hat / sqrt(ssq_dndgnX_beta_hat))))
    
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      match_on_more = match_on_more,
      regress_on_X = regress_on_X,
      infrence = "discordant no diffs",
      beta_hat_T = dndgnX_beta_hat,
      ssq_hat_T = ssq_dndgnX_beta_hat,
      pval = dndgnX_pval
    ))
    
    diffs_y_0_1 = as.numeric(diffs_y == 1)
    discordant_glm_no_X = summary(glm(diffs_y_0_1 ~ 0 + diffs_treatment, family = "binomial"))$coefficients[1,c(1,2)]
    dgnX_beta_hat = discordant_glm_no_X[1]
    ssq_dgnX_beta_hat = discordant_glm_no_X[2]^2
    dgnX_pval = 2 * pnorm(min(c(-1,1) * (dgnX_beta_hat / sqrt(ssq_dgnX_beta_hat))))
    
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      match_on_more = match_on_more,
      regress_on_X = regress_on_X,
      infrence = "conditional logit",
      beta_hat_T = dgnX_beta_hat,
      ssq_hat_T = ssq_dgnX_beta_hat,
      pval = dgnX_pval
    ))
    
    logit_model_glm_no_X = summary(glm(y ~ w, family = "binomial"))$coefficients[2,c(1,2)]
    lgnX_beta_hat = logit_model_glm_no_X[1]
    ssq_lgnX_beta_hat = logit_model_glm_no_X[2]^2
    lgnX_pval = 2 * pnorm(min(c(-1,1) * (lgnX_beta_hat / sqrt(ssq_lgnX_beta_hat))))
    
    res = rbind(res, data.frame(
      n = n,
      beta_T = beta_T,
      match_on_more = match_on_more,
      regress_on_X = regress_on_X,
      infrence = "logit",
      beta_hat_T = lgnX_beta_hat,
      ssq_hat_T = ssq_lgnX_beta_hat,
      pval = lgnX_pval
    ))
  }
  
  if (FALSE) {
    # discordant_no_diffs_model_glm = summary(glm(y[discordant_idx] ~ w[discordant_idx] + X[discordant_idx], family = "binomial"))$coefficients[2,c(1,2)]
    # dndg_beta_hat = discordant_no_diffs_model_glm[1]
    # ssq_dndg_beta_hat = discordant_no_diffs_model_glm[2]^2
    # dndg_pval = 2 * pnorm(min(c(-1,1) * (dndg_beta_hat / sqrt(ssq_dndg_beta_hat))))
    # 
    # res = rbind(res, data.frame(
    #   n = n,
    #   beta_T = beta_T,
    #   infrence = "discordant no diffs glm",
    #   beta_hat_T = dndg_beta_hat,
    #   ssq_hat_T = ssq_dndg_beta_hat,
    #   pval = dndg_pval
    # ))
    
    # discordant_glm = summary(glm(diffs_y_0_1 ~ 0 + diffs_treatment + diffs_X, family = "binomial"))$coefficients[1,c(1,2)]
    # dg_beta_hat = discordant_glm[1]
    # ssq_dg_beta_hat = discordant_glm[2]^2
    # dg_pval = 2 * pnorm(min(c(-1,1) * (dg_beta_hat / sqrt(ssq_dg_beta_hat))))
    # 
    # res = rbind(res, data.frame(
    #   n = n,
    #   beta_T = beta_T,
    #   match_on_more = match_on_more,
    #   infrence = "discordant glm",
    #   beta_hat_T = dg_beta_hat,
    #   ssq_hat_T = ssq_dg_beta_hat,
    #   pval = dg_pval
    # ))
    
    # m_beta_hat = mixed_model$mixed_betaT
    # ssq_m_beta_hat = mixed_model$mixed_ssq_b
    # m_pval = 2 * pnorm(min(c(-1,1) * (m_beta_hat / sqrt(ssq_m_beta_hat))))
    # 
    # res = rbind(res, data.frame(
    #   n = n,
    #   beta_T = beta_T,
    #   match_on_more = match_on_more,
    #   infrence = "mixed",
    #   beta_hat_T = m_beta_hat,
    #   ssq_hat_T = ssq_m_beta_hat,
    #   pval = m_pval
    # ))
    
    # c_beta_hat = mixed_model$concordnat_betaT
    # ssq_c_beta_hat = mixed_model$concordnat_ssq_b
    # c_pval = 2 * pnorm(min(c(-1,1) * (c_beta_hat / sqrt(ssq_c_beta_hat))))
    # 
    # res = rbind(res, data.frame(
    #   n = n,
    #   beta_T = beta_T,
    #   infrence = "concordant",
    #   beta_hat_T = c_beta_hat,
    #   ssq_hat_T = ssq_c_beta_hat,
    #   pval = c_pval
    # ))
    
    # logit_model_glm = summary(glm(y ~ w + X, family = "binomial"))$coefficients[2,c(1,2)]
    # lg_beta_hat = logit_model_glm[1]
    # ssq_lg_beta_hat = logit_model_glm[2]^2
    # lg_pval = 2 * pnorm(min(c(-1,1) * (lg_beta_hat / sqrt(ssq_lg_beta_hat))))
    # 
    # res = rbind(res, data.frame(
    #   n = n,
    #   beta_T = beta_T,
    #   infrence = "logit glm",
    #   beta_hat_T = lg_beta_hat,
    #   ssq_hat_T = ssq_lg_beta_hat,
    #   pval = lg_pval
    # ))
  }
  return(res)
}


# for (j in 1:1000) {
#   print(j)
#   row = params[j,]
#   print(row)
# 
#   i = row$i
#   beta_T = row$beta_T
#   n = row$n
#   match_on_more = row$match_on_more
#   regress_on_X = row$regress_on_X
# 
#   tryCatch({
#     out = run_simulation(i, beta_T, n, match_on_more, regress_on_X)
#   }, error = function(e) {
#     cat(glue::glue("Error in sim i={i}: {e$message}"), '\n')
#     out = NULL
#   })
#   #print(out)
# }


handlers(global = TRUE)
handlers("txtprogressbar")

registerDoFuture()
plan(multisession, workers = num_cores)


start_time = Sys.time()

for (e_nsim in 1:external_nsim){
  handlers(global = TRUE)
  handlers("txtprogressbar")
  
  registerDoFuture()
  plan(multisession, workers = num_cores)
  
  tryCatch({
    with_progress({
      prog = progressor(along = 1:nrow(params))
      results = foreach(row = iter(params, by = "row"),
                        .combine = rbind,
                        .packages = c("clogitR", "nbpMatching", "data.table",
                                      "dplyr", "MASS", "Rcpp", "rstanarm")) %dorng% {

        #set.seed(row$i + (as.numeric(Sys.time())))
        i = row$i
        beta_T = row$beta_T
        n = row$n
        match_on_more = row$match_on_more
        regress_on_X = row$regress_on_X
        #cat(glue::glue("Running i={i}"), '\n')
        res = tryCatch({
          out = run_simulation(i, beta_T, n, match_on_more, regress_on_X)
          #cat("Successfully ran simulation")
          prog()
          out
        }, error = function(e) {
          cat(glue::glue("Error in sim i={i}: {e$message}"), '\n')
          prog()  # still update progress bar even if it fails
          NULL    # return NULL if failed, will be dropped in rbind
        })
      }
    })
    write.csv(results, file = paste0("C:/temp/clogitR_kap_test_3/100_", e_nsim, ".csv"), row.names = FALSE)
    rm(results); gc()
  }, error = function(e) {
    cat(glue::glue("Error in e_nsim ={e_nsim}: {e$message}"), '\n')
  })
  plan(sequential)
}

# for (e_nsim in 1:external_nsim) {
# 
#   with_progress({
#     prog = progressor(along = 1:nrow(params))
# 
#     results = foreach(row = iter(params, by = "row"),
#                       .combine = rbind,
#                       .packages = c("clogitR", "nbpMatching", "data.table",
#                                     "dplyr", "MASS", "Rcpp", "rstanarm")) %dorng% {
# 
#       # extract parameters
#       i = row$i
#       beta_T = row$beta_T
#       n = row$n
#       match_on_more = row$match_on_more
#       regress_on_X = row$regress_on_X
# 
#       # retry mechanism for random Stan failures
#       attempts = 0
#       max_attempts = 3
#       out = NULL
# 
#       while (attempts < max_attempts && is.null(out)) {
#         attempts = attempts + 1
# 
#         out = tryCatch({
#           # prevent nested Stan parallelization
#           #options(mc.cores = 1)
# 
#           # timeout in case Stan hangs
#           R.utils::withTimeout({
#             run_simulation(i, beta_T, n, match_on_more, regress_on_X)
#           }, timeout = 180, onTimeout = "silent")  # adjust seconds as needed
# 
#         }, error = function(e) {
#           message(glue("Error in i={i}, attempt {attempts}: {e$message}"))
#           NULL
#         })
# 
#         if (is.null(out)) Sys.sleep(runif(1, 1, 3))  # random small pause before retry
#       }
# 
#       prog()  # update progress bar
#       out  # return result (NULLs are dropped automatically by rbind)
#     }
#   })
# 
#   write.csv(results, file = paste0("C:/temp/clogitR_kap_test/500_", e_nsim, ".csv"), row.names = FALSE)
#   rm(results); gc()
# 
# }


plan(sequential)
end_time = Sys.time()


results = read.csv("C:/temp/clogitR_kap_test_3/100_1.csv")

sum = 0
for (i in 2:475) {
  file_path <- paste0("C:/temp/clogitR_kap_test_3/100_", i, ".csv")
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
  results = rbind(results, read.csv(paste0("C:/temp/clogitR_kap_test/100_", i, ".csv")))
}
500 * i
results$X = NULL

res_mod = results %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(beta_T, match_on_more, infrence, n, regress_on_X) %>%
  summarize(
    num_na = sum(is.na(pval)), 
    mse = mean(sq_err, na.rm = TRUE),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    mean_betaT_hat = mean(beta_hat_T, na.rm = TRUE),
    mean_sq_betaT_hat = mean(sqrt(ssq_hat_T), trim = 0.001, na.rm = TRUE),
    .groups = "drop")


write.csv(res_mod, file = "C:/temp/clogitR_kap_test_3/combined_6000.csv", row.names = FALSE)

  
  
  
Sigma = 1 * (matrix(0.75, nrow = 6, ncol = 6) + diag(1 - 0.75, 6))
X = data.table(MASS::mvrnorm(n, rep(1, 6), Sigma))

df = data.frame(cbind(id = 1:n, X))
df.dist = gendistance(df[, -1], idcol = 1)
df.mdm = distancematrix(df.dist)
df.match = nonbimatch(df.mdm)

T_inx = df.match$halves[,2]
w = as.numeric(1:n %in% T_inx)

strat = rep(0, n)

for(pair_idx in seq_len(nrow(df.match$halves))) {
  strat[df.match$halves[pair_idx, 2]] = pair_idx
  strat[df.match$halves[pair_idx, 4]] = pair_idx
}
beta_0 = -1
beta_x = c(1, 2, -2, -3, 3, 0)
probs = 1 / (1 + exp(-(beta_0 + (as.matrix(X) %*% beta_x) + beta_T * w)))
y = rbinom(n, 1, probs)

df = data.frame(
  p = probs,
  y = y
)

ggplot(df, aes(x = p, fill = factor(y))) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
  labs(x = "Predicted probability", fill = "Outcome") +
  scale_fill_manual(values = c("0" = "red", "1" = "blue")) +
  theme_minimal()

#drop the X to one col:
X = X[,1]
  
  
  
pacman::p_load(fastLogisticRegressionWrap)

n = 100
p = 4
alpha = 0.05

x = runif(n * p, -1, 1)
if (p == 1){
  x = sort(x)
}
X = cbind(1, matrix(x, ncol = p))
beta_x_vec = rep(3, p)
beta_0 = -1
beta_0_plus_beta_x_times_x = X %*% c(beta_0, beta_x_vec)
beta_T = 0.75
logit = function(x){1 / (1 + exp(-x))}

Nsim = 10000
num_disc = array(NA, Nsim)
pval_x = matrix(NA, Nsim, p)
pval_w_x = array(NA, Nsim)
pval_w = array(NA, Nsim)
for (nsim in 1 : Nsim){
  w = if (p == 1){
    rep(sample(c(0, 1)), times = n / 2)
  } else {
    sample(c(rep(0, n / 2), rep(1, n / 2)))
  }
  
  prob_y_eq_one = logit(beta_0_plus_beta_x_times_x + beta_T * w)
  y = rbinom(n, 1, prob_y_eq_one)
  
  obj = fast_logistic_regression(X, y,           do_inference_on_var = "all")
  pval_x[nsim, ] = obj$approx_pval[-1]
  obj = fast_logistic_regression(cbind(X, w), y, do_inference_on_var = p + 2)
  pval_w_x[nsim] = obj$approx_pval[p + 2]
  obj = fast_logistic_regression(cbind(1, w), y, do_inference_on_var = 2)
  pval_w[nsim] = obj$approx_pval[2]
  # y_and_w = matrix(y, nrow = n/2, byrow = TRUE)
  # num_disc[nsim] = sum(rowSums(y_and_w) == 1)
}
# cbind(y_and_w, matrix(p, nrow = n/2, byrow = TRUE))
# summary(pval_x)
# summary(pval_w_x)
# summary(pval_w)
colMeans(pval_x < alpha)
mean(pval_w_x[!is.nan(pval_w_x)] < alpha)
mean(pval_w < alpha)

df = data.frame(
  p = prob_y_eq_one,
  y = y
)

ggplot(df, aes(x = p, fill = factor(y))) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
  labs(x = "Predicted probability", fill = "Outcome") +
  scale_fill_manual(values = c("0" = "red", "1" = "blue")) +
  theme_minimal()
  




