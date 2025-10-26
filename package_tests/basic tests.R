#pacman::p_load(clogitR, dplyr, data.table, ggplot2, nbpMatching, survival)
devtools::load_all("C:/Users/Jacob/clogitR_package_and_paper_repo/clogitR")
rm(list = ls())
options(error = recover)

#n = 50
#p = 2
mu_x = 1
sigma_x = 1
sigma_e = 1
#nsim_exact_test = 501
num_cores = 1
Nsim = 10

#build mvnp covariates
set.seed(1)


all_betas_and_correlations = list()
all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 2, rho = 0,     betas = c(1, -3, 1.5, 0, 0))))
all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 2, rho = 0.75,  betas = c(1, -3, 1.5, 0, 0))))
all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 2, rho = 0,     betas = c(0.5, -2, 1, 0, 0))))
all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 2, rho = 0.75,  betas = c(0.5, -2, 1, 0, 0))))
all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 4, rho = 0,     betas = c(1, -2, 1, 1.5, -2))))
all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 4, rho = 0.75,  betas = c(1, -2, 1, 1.5, -2))))
# all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 20, rho = 0,     betas = c(1, -2, 1, 1.5, -2,
#                                                                                                     0, 0, -1, -3, 2,
#                                                                                                     1, -2, 0, 1.4, 0,
#                                                                                                     -1.3, -2, 3, 1, 1
# ))))
# all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 20, rho = 0.75,  betas = c(1, -2, 1, 1.5, -2,
#                                                                                                     0, 0, -1, -3, 2,
#                                                                                                     1, -2, 0, 1.4, 0,
#                                                                                                     -1.3, -2, 3, 1, 1
# ))))

# plots_list = list()
# plot_index = 1

comparison = data.frame(
  n = numeric(),
  beta_T = numeric(),
  p = numeric(),
  num_discord = numeric(),
  Rcpp_beta_hat_T = numeric(),
  Rcpp_sse_beta_T = numeric(),
  clogit_beta_hat_T = numeric(),
  clogit_sse_beta_T = numeric(),
  glm_beta_hat_T = numeric(),
  glm_sse_beta_T = numeric(),
  swaped = numeric()
)
to_few_discordant = list()

res = data.frame(
  betas = character(),
  beta_T = numeric(),
  n = numeric(),
  rho = numeric(),
  infrence = character(),
  beta_hat_T = numeric(),
  pval = numeric()
)
times = rep(NA, Nsim)
for (nsim in 1 : Nsim){
  start_time = Sys.time()
  cat ("nsim:", nsim, "/", Nsim, "\n")
  
  tryCatch({
  for (beta_T in c(0,1)){
    for (n in c(500, 1000)){ #50, 100, 500, 1000
      #beta_T = 0; n = 500; all_betas_and_correlation = all_betas_and_correlations[[1]]
      errors = rnorm(n, 0, sigma_e)
      for (all_betas_and_correlation in all_betas_and_correlations){
        probs = rep(NA, n)
        y = rep(NA, n)
        p = all_betas_and_correlation[["p"]]
        betas = all_betas_and_correlation[["betas"]]
        rho = all_betas_and_correlation[["rho"]] 
        
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
        
        for(i in seq_len(nrow(df.match$halves))) {
          strat[df.match$halves[i,2]] = i
          strat[df.match$halves[i,4]] = i
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
        
        
        # df = data.frame(
        #   p = probs,
        #   y = y
        # )
        # plots_list[[plot_index]] =
        #   ggplot(df, aes(x = p, fill = factor(y))) +
        #     geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
        #     labs(x = "Predicted probability", fill = "Outcome") +
        #     scale_fill_manual(values = c("0" = "red", "1" = "blue")) +
        #     theme_minimal()
        # plot_index = plot_index + 1
        
        #test out model, clogit from survival, and simple glm
        
        #clogitR_model = clogitR(response = y, data = X, treatment = w_t, strata = strat, use_concordant_pairs = TRUE, try_to_stableize = TRUE)
        #clogit(y ~ w_t + . + strata(strat), data = X)
        #clogitR_model$coefficients_and_variances
        #summary(clogitR_model)
        
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
        
      }
    }
  }
  times[nsim] = Sys.time() - start_time
  print(Sys.time() - start_time)
  }, error = function(e) {
    message(paste("⚠️ Error in nsim =", nsim, ":", conditionMessage(e)))
    times[nsim] = NA
  })
}
res_mod = res %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(betas, rho, beta_T, n, infrence) %>%
  summarize(
    num_na = sum(is.na(pval)), 
    mse = mean(sq_err, na.rm = TRUE),
    percent_reject = sum(rej, na.rm = TRUE) / (n() - num_na),
    .groups = "drop")


pacman::p_load(microbenchmark)
Rcpp::sourceCpp("C:/Users/Jacob/matching_on_the_fly_designs_R_package_and_paper_repr/package_tests/fast_conditional_logistic_regr.cpp")

df = matched_inference$get_profiling_data()$clogit_df

X_discord = matched_inference$get_profiling_data()$X_matched_diffs_discordant
y_discord = matched_inference$get_profiling_data()$y_matched_diffs_discordant


microbenchmark(
  clogit = clogit(y ~ w + .-str + strata(str), data = df),
  Rcpp = fast_conditional_logistic_regression_with_var_cpp(cbind(1, X_discord), y_discord, j = 1),
  glm = glm(y_discord_1_0 ~ X_discord, family = "binomial"),
  times = 5000L
)

clogit_mod = clogit(y ~ w + .-str + strata(str), data = df)
Rcpp_mod = fast_conditional_logistic_regression_with_var_cpp(cbind(1, X_discord), y_discord, j = 1)

y_discord_1_0 = as.numeric(y_discord == 1)
glm_mod = glm(y_discord_1_0 ~ X_discord, family = "binomial")






comparison$diff_beta_hat_t = comparison$Rcpp_beta_hat_T - comparison$clogit_beta_hat_T
comparison$diff_sse_hat_t = comparison$Rcpp_sse_beta_T - comparison$clogit_sse_beta_T

ggplot(data = comparison, aes(x = Rcpp_beta_hat_T, y = clogit_beta_hat_T, colour = factor(beta_T), shape = factor(swaped))) +
  geom_point() +
  facet_wrap(~ n)

ggplot(data = comparison, aes(x = Rcpp_beta_hat_T, y = glm_beta_hat_T, colour = factor(beta_T), shape = factor(swaped))) +
  geom_point() +
  facet_wrap(~ n)

ggplot(data = comparison, aes(x = num_discord, fill = factor(swaped))) + 
  geom_histogram() +
  facet_wrap(~ p)


ggplot(data = comparison, aes(x = Rcpp_sse_beta_T, y = clogit_sse_beta_T, colour = factor(beta_T), shape = factor(swaped))) +
  geom_point() +
  facet_wrap(~ n)

ggplot(data = comparison, aes(x = Rcpp_sse_beta_T, y = clogit_sse_beta_T, colour = factor(beta_T), shape = factor(swaped))) +
  geom_point() +
  xlim(0, 1e1) +
  ylim(0, 1e1)

















