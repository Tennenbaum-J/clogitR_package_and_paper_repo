pacman::p_load(clogitR, dplyr, data.table, ggplot2, nbpMatching)
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
all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 20, rho = 0,     betas = c(1, -2, 1, 1.5, -2,
                                                                                                    0, 0, -1, -3, 2,
                                                                                                    1, -2, 0, 1.4, 0,
                                                                                                    -1.3, -2, 3, 1, 1
))))
all_betas_and_correlations = c(all_betas_and_correlations, list(list(p = 20, rho = 0.75,  betas = c(1, -2, 1, 1.5, -2,
                                                                                                    0, 0, -1, -3, 2,
                                                                                                    1, -2, 0, 1.4, 0,
                                                                                                    -1.3, -2, 3, 1, 1
))))

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
  design = character(),
  infrence = character(),
  beta_hat_T = numeric(),
  pval = numeric()
)
times = rep(NA, Nsim)
for (nsim in 1 : Nsim){
  start_time = Sys.time()
  cat ("nsim:", nsim, "/", Nsim, "\n")
  
  #tryCatch({
  for (beta_T in c(0,1)){
    for (n in c(50, 100, 500, 1000)){ #50, 100, 500, 1000
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
        
        clogitR_model = clogitR(y ~ treatment(w_t) + X + strata(strat), use_concordant_pairs = TRUE, try_to_stableize = TRUE)
        
        
        
        matched_inference = SeqDesignInferenceIncidConditionalLogRegrKK$new(seq_des_obj, num_cores = 1, convex_flag = TRUE, verbose = FALSE)
        mixed_inference = SeqDesignInferenceIncidConditionalLogRegrKK$new(seq_des_obj, num_cores = 1, convex_flag = FALSE, verbose = FALSE)
        
        #cat("bata  : ", sprintf("%.3f", beta_T), ' ', sprintf("%.3f", betas[1]), ' ', sprintf("%.3f", betas[2]), '\n')
        matched_beta_hat_T = matched_inference$compute_treatment_estimate()
        matched_pval = matched_inference$compute_mle_two_sided_pval_for_treatment_effect()
        
        mixed_beta_hat_T = mixed_inference$compute_treatment_estimate()
        mixed_pval = mixed_inference$compute_mle_two_sided_pval_for_treatment_effect()
        #cat('\n')
        
        if(!is.null(matched_inference$get_graphing_data()$Rcpp_beta_hat_T)){
          comparison = rbind(comparison, data.frame(
            n = n,
            beta_T = beta_T,
            p = p,
            num_discord = matched_inference$get_graphing_data()$num_discord,
            Rcpp_beta_hat_T = matched_inference$get_graphing_data()$Rcpp_beta_hat_T,
            Rcpp_sse_beta_T = matched_inference$get_graphing_data()$Rcpp_sse_beta_T,
            clogit_beta_hat_T = matched_inference$get_graphing_data()$clogit_beta_hat_T,
            clogit_sse_beta_T = matched_inference$get_graphing_data()$clogit_sse_beta_T,
            glm_beta_hat_T = matched_inference$get_graphing_data()$glm_beta_hat_T,
            glm_sse_beta_T = matched_inference$get_graphing_data()$glm_sse_beta_T,
            swaped = matched_inference$get_graphing_data()$swaped
          ))
          to_few_discordant = c(to_few_discordant, 0)
        } else {
          to_few_discordant = c(to_few_discordant, 1)
        }
        res = rbind(res, data.frame(
          betas = paste0(betas, collapse=""),
          beta_T = beta_T,
          n = n,
          rho = rho,
          design = d,
          infrence = "matched",
          beta_hat_T = matched_beta_hat_T,
          pval = matched_pval
        ))
        res = rbind(res, data.frame(
          betas = paste0(betas, collapse=""),
          beta_T = beta_T,
          n = n,
          rho = rho,
          design = d,
          infrence = "mixed",
          beta_hat_T = mixed_beta_hat_T,
          pval = mixed_pval
        ))
        
      }
    }
  }
  times[nsim] = Sys.time() - start_time
  print(Sys.time() - start_time)
  # }, error = function(e) {
  #   message(paste("⚠️ Error in nsim =", nsim, ":", conditionMessage(e)))
  #   times[nsim] = NA  
  # })
}
res_mod = res %>%
  mutate(sq_err = (beta_hat_T - beta_T)^2, rej = pval < 0.05) %>%
  group_by(betas, rho, beta_T, n, design, infrence) %>%
  summarize(mse = mean(sq_err), percent_reject = sum(rej) / n())


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


