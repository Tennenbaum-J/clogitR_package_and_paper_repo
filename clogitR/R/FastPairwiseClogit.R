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
    discordant_model = fast_conditional_logistic_regression_with_var_cpp(cbind(discordnat_Treatmentdiffs, discordant_Xdiffs),
                                                                         discordant_ydiffs, j = 1)
  }
  if (!is.null(concordnat_X) & !is.null(concordnat_y) & !is.null(concordnat_Treatment)) {
    concordnat_y_1_n1 = 2 * concordnat_y - 1
    concordant_model = fast_conditional_logistic_regression_with_var_cpp(cbind(1, concordnat_Treatment, concordnat_X),
                                                                         concordnat_y_1_n1, j = 2)
  }
  ret = list()
  if (!is.null(discordant_model)){
    ret$discordnat_betaT = discordant_model$b[1]
    ret$discordnat_ssq_b = discordant_model$ssq_b_j
  }
  if (!is.null(concordant_model)){
    ret$concordnat_betaT = concordant_model$b[2]
    ret$concordnat_ssq_b = concordant_model$ssq_b_j
  }
  if (!is.null(discordant_model) && !is.null(concordant_model)){
    w_star = concordant_model$ssq_b_j / (concordant_model$ssq_b_j + discordant_model$ssq_b_j)
    ret$mixed_betaT = discordant_model$b[1] * w_star + concordant_model$b[2] * (1 - w_star)
    ret$mixed_ssq_b = discordant_model$ssq_b_j * w_star
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
    concordnat_y_1_n1 = 2 * concordnat_y - 1
    concordant_model = fast_conditional_logistic_regression_with_var_cpp(cbind(1, concordnat_Treatment, concordnat_X),
                                                                         concordnat_y_1_n1, j = -1)
  }
  
  
  if (!is.null(discordant_Xdiffs) && !is.null(discordant_ydiffs) && !is.null(discordnat_Treatmentdiffs)) {
    # Check if we have a concordant model to use as prior
    if (!is.null(concordant_model)) {
      # Extract prior from concordant model (excluding intercept)
      prior_cov_no_intercept = c(20, sqrt(concordant_model$ssq[-(1:2)]))
      prior_means = c(0, concordant_model$b[-(1:2)])
      
      # Convert to 0/1 if needed
      discordant_ydiffs_0_1 = ifelse(discordant_ydiffs == -1, 0, 1)
      
      # Use rstanarm with informative prior
      discordant_model = rstanarm::stan_glm(
        discordant_ydiffs_0_1 ~ 0 + discordnat_Treatmentdiffs + discordant_Xdiffs,
        family = binomial(link = "logit"),
        prior = rstanarm::normal(location = prior_means, scale = prior_cov_no_intercept),
        data = data.frame(discordant_ydiffs_0_1, discordnat_Treatmentdiffs, discordant_Xdiffs),
        refresh = 0  # Suppress Stan output
      )
    } else {
      warning("No concordant model available for prior. Skipping discordant model.")
    }
  }
  
  
  # #the concordant model has an intersect. We don't use that in the discordnat model, so remove it from the prior
  # prior_cov_no_intercept = sqrt(concordant_model$ssq[-1])
  # 
  # # Get prior means (excluding intercept)
  # prior_means = concordant_model$b[-1]
  # 
  # discordant_ydiffs_0_1 = as.numeric(discordant_ydiffs == 1)
  # discordant_model = rstanarm::stan_glm(
  #   discordant_ydiffs_0_1 ~ 0 + discordnat_Treatmentdiffs + discordant_Xdiffs,
  #   family = binomial(link = "logit"),
  #   prior = rstanarm::normal(location = prior_means, scale = prior_cov_no_intercept),
  #   data = data.frame(discordant_ydiffs_0_1, discordnat_Treatmentdiffs, discordant_Xdiffs), refresh=0
  # )
  # #, refresh=0
  
  ret = list()
  if (!is.null(discordant_model)){
    ret$discordnat_betaT = unname(coef(discordant_model)[1])
    posterior_draws = as.matrix(discordant_model)
    ret$discordnat_ssq_b = var(posterior_draws[, 1])
  }
  
  if (!is.null(concordant_model)){
    ret$concordnat_betaT = concordant_model$b[2]
    ret$concordnat_ssq_b = concordant_model$ssq[2]
  }
  
  if (!is.null(discordant_model) && !is.null(concordant_model)){
    w_star = ret$concordnat_ssq_b / (ret$concordnat_ssq_b + ret$discordnat_ssq_b)
    ret$mixed_betaT = ret$discordnat_betaT * w_star + ret$concordnat_betaT * (1 - w_star)
    ret$mixed_ssq_b = ret$discordnat_ssq_b * w_star
  }
  
  return(ret)
}