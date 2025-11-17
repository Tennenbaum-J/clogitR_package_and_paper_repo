#' Fast Conditional logistic regression with prepared data without any checks
#'
#' @description
#' Computes the treatment effect and its squared standard error for the discordant and concordant pairs.
#' It then returns a list containing theses values as well as the mixed estimates.
#' 
#' This is the fast, very bare-bones, version of the model, primary meant for power simulations or other situations where
#' the overhead of object creation is too cumbersome. As such there are no checks, validations, or warnings, in this model.
#' It will break if there are not enough renovates, or the length of vectors don't match, and doesn't make any attempt to 
#' stabilize potentially unstable results. USE AT YOUR OWN RISK.
#' 
#' For the standard model, see `clogitR`.
#' 
#' 
#' @param discordant_Xdiffs A fully prepared data matrix of the differences of the discordant matched pairs.
#' @param discordant_ydiffs A fully prepared binary (-1,1) vector of the differences in response the discordant matched pairs.
#' @param discordnat_Treatmentdiffs A fully prepared binary (-1,1) vector of the differences in treatment the discordant matched pairs.
#' @param concordnat_X A fully prepared data matrix of the concordant matched pairs or reservoir.
#' @param concordnat_y A fully prepared binary (0,1) vector of the response the concordant matched pairs or reservoir.
#' @param concordnat_Treatment A fully prepared binary (0,1) vector of the treatment the concordant matched pairs or reservoir.
#' @return A list containing the estimated treatment effect and squared standard error for each of the applicable models (discordant, concordant, mixed).
#' @export
fastClogit = function(discordant_Xdiffs = NULL,
                      discordant_ydiffs = NULL,
                      discordnat_Treatmentdiffs = NULL,
                      concordnat_X = NULL,
                      concordnat_y = NULL,
                      concordnat_Treatment = NULL) {

  discordant_model = NULL
  concordant_model = NULL
  
  if (!is.null(discordant_Xdiffs) & !is.null(discordant_ydiffs) & !is.null(discordnat_Treatmentdiffs)) {
    #discordant_model = fast_conditional_logistic_regression_with_var_cpp(cbind(discordnat_Treatmentdiffs, discordant_Xdiffs),
    #                                                                     discordant_ydiffs, j = 1)
    discordant_ydiffs_0_1 = ifelse(discordant_ydiffs == -1, 0, 1)
    discordant_model = summary(glm(discordant_ydiffs_0_1 ~ 0 + discordnat_Treatmentdiffs + discordant_Xdiffs, family = "binomial"))$coefficients[1,c(1,2)]
  }
  if (!is.null(concordnat_X) & !is.null(concordnat_y) & !is.null(concordnat_Treatment)) {
    #concordnat_y_1_n1 = 2 * concordnat_y - 1
    #concordant_model = fast_conditional_logistic_regression_with_var_cpp(cbind(1, concordnat_Treatment, concordnat_X),
    #                                                                     concordnat_y_1_n1, j = 2)
    concordant_model = summary(glm(concordnat_y ~ concordnat_Treatment + concordnat_X, family = "binomial"))$coefficients[2,c(1,2)]
  }
  ret = list()
  if (!is.null(discordant_model)){
    # ret$discordnat_betaT = discordant_model$b[1]
    # ret$discordnat_ssq_b = discordant_model$ssq_b_j
    ret$discordnat_betaT = discordant_model[1]
    ret$discordnat_ssq_b = discordant_model[2]^2
  }
  if (!is.null(concordant_model)){
    # ret$concordnat_betaT = concordant_model$b[2]
    # ret$concordnat_ssq_b = concordant_model$ssq_b_j
    ret$concordnat_betaT = concordant_model[1]
    ret$concordnat_ssq_b = concordant_model[2]^2
  }
  if (!is.null(discordant_model) && !is.null(concordant_model)){
    w_star = concordant_model[2]^2 / (concordant_model[2]^2 + discordant_model[2]^2)
    ret$mixed_betaT = discordant_model[1] * w_star + concordant_model[1] * (1 - w_star)
    ret$mixed_ssq_b = discordant_model[2]^2 * w_star
  }
  return(ret)
}







#' @export
bayesianClogit = function(discordant_Xdiffs = NULL,
                          discordant_ydiffs = NULL,
                          discordnat_Treatmentdiffs = NULL,
                          concordnat_X = NULL,
                          concordnat_y = NULL,
                          concordnat_Treatment = NULL) {
  
  discordant_model = NULL
  concordant_model = NULL
  
  if (!is.null(concordnat_X) && !is.null(concordnat_y) && !is.null(concordnat_Treatment)) {
    if (any(is.na(concordnat_X)) || any(is.na(concordnat_y)) || any(is.na(concordnat_Treatment))) {
      warning("NA/NaN values detected in concordant data. Skipping concordant model.")
    } else {
      concordnat_y_1_n1 = 2 * concordnat_y - 1
      concordant_model = tryCatch({
        summary(glm(concordnat_y ~ concordnat_Treatment + concordnat_X, family = "binomial"))$coefficients
        #summary(lm(concordnat_y ~ concordnat_Treatment + concordnat_X))$coefficients
        #fast_conditional_logistic_regression_with_var_cpp(cbind(1, concordnat_Treatment, concordnat_X),
        #                                                  concordnat_y_1_n1, j = -1)
      }, error = function(e) {
        warning(sprintf("fast_conditional_logistic_regression_with_var_cpp failed: %s", e$message))
        NULL
      })
    }
  }
  
  if (!is.null(discordant_Xdiffs) && !is.null(discordant_ydiffs) && !is.null(discordnat_Treatmentdiffs)) {
    if (any(is.na(discordant_Xdiffs)) || any(is.na(discordant_ydiffs)) || any(is.na(discordnat_Treatmentdiffs))) {
      warning("NA/NaN values detected in discordant data. Skipping discordant model.")
    } else if (!is.null(concordant_model)) {
      prior_means = concordant_model[-1,1]
      prior_cov_no_intercept = min(concordant_model[-1,2], 20)
      if (sum(prior_cov_no_intercept) == 20 * length(prior_cov_no_intercept)) {
        # discordant_ydiffs_0_1 = ifelse(discordant_ydiffs == -1, 0, 1)
        # discordant_model = summary(glm(discordant_ydiffs_0_1 ~ 0 + discordnat_Treatmentdiffs + discordant_Xdiffs, family = "binomial"))$coefficients
        # ret = list(
        #   discordnat_betaT = discordant_model[1,1],
        #   discordnat_ssq_b = discordant_model[1,2]
        # )
        return(NULL)
      }
      
      if (any(is.na(prior_cov_no_intercept)) || any(is.na(prior_means)) || 
          any(!is.finite(prior_cov_no_intercept)) || any(!is.finite(prior_means))) {
        warning("Prior contains NA/NaN/Inf. Skipping discordant model.")
      } else {
        discordant_ydiffs_0_1 = ifelse(discordant_ydiffs == -1, 0, 1)
        
        discordant_model = tryCatch({
          rstanarm::stan_glm(
            discordant_ydiffs_0_1 ~ 0 + discordnat_Treatmentdiffs + discordant_Xdiffs,
            family = binomial(link = "logit"),
            prior = rstanarm::normal(location = prior_means, scale = prior_cov_no_intercept),
            data = data.frame(discordant_ydiffs_0_1, discordnat_Treatmentdiffs, discordant_Xdiffs),
            refresh = 0
          )
        }, error = function(e) {
          warning(sprintf("stan_glm failed: %s", e$message))
          NULL
        })
      }
    } else {
      warning("No valid concordant model available for prior. Skipping discordant model.")
    }
  }
  
  ret = list()
  if (!is.null(discordant_model)){
    ret$discordnat_betaT = unname(coef(discordant_model)[1])
    posterior_draws = as.matrix(discordant_model)
    ret$discordnat_ssq_b = var(posterior_draws[, 1])
  }
  
  return(ret)
}