#' Initialize a new clogitR model
#'
#' @param response A binary (0,1) vector containing the response of each subject.
#' @param data A data.frame, data.table, or model.matrix containing the variables.
#' @param treatment A binary (0,1) vector containing weather each subject is treatment or control.
#' @param strata Numeric vector indexing the matched pairs and the reservoir. Entrees indexed with zero, are added to the reservoir. Entrees indexed with the same number are part of the same strata.
#' @param use_concordant_pairs Flag indicating whether to include the concordant pairs and reservoir in the model.
#' @param concordant_method The method to use for fitting the concordant pairs and reservoir. Options are "GLM", "GEE", and "GLMM".
#' @param try_to_stableize Flag indicating whether the model should decompose the discordant pairs to increase stability if the model is unstable.
#' @return A list of class `clogitR`.
#' @export
clogitR = function(response = NULL,
                   data = NULL,
                   treatment = NULL,
                   strata = NULL,
                   use_concordant_pairs = TRUE,
                   concordant_method = "GEE",
                   try_to_stableize = TRUE) {
  
  # ------------------------------------------------------------------------
  # 1. Input Validation and Pre-processing
  # ------------------------------------------------------------------------
  if (test_data_frame(data, types = c("numeric", "integer", "factor", "logical"))) {
    data_mat = model.matrix(~0+., data = data)
  } else if (test_matrix(data, mode = "numeric")) {
    if (sd(data[, 1]) == 0) { data[, 1] = NULL }
    data_mat = as.matrix(data)
  } else {
    assert(
      check_data_frame(data),
      check_matrix(data, mode = "numeric")
    )
  }
  
  n = nrow(data_mat)
  data_name = deparse(substitute(data))
  
  assertFlag(use_concordant_pairs)
  assertString(concordant_method, c("GLM", "GEE", "GLMM"))
  assertFlag(try_to_stableize)

  
  
  assertNumeric(response, lower = 0, upper = 1, any.missing = FALSE, len = n)
  response_name = deparse(substitute(response))
  
  treatment_name = NULL
  if (!is.null(treatment)) {
    assertNumeric(treatment, lower = 0, upper = 1, any.missing = FALSE, len = n)
    if (!all(treatment %in% c(0,1))) { stop("Treatment must be binary 0 or 1.") }
    treatment_name = deparse(substitute(treatment))
  }
  
  assertNumeric(strata, any.missing = FALSE, len = n)
  strata_name = deparse(substitute(strata))
  
  # Capture terms for summary/printing usage if possible (though we used matrix for fitting)
  # If data was a dataframe, we can't easily reconstruction 'terms' from the matrix unless we kept the df.
  # But the original R6 used `terms(response ~ ., data = self$data)` which seems to imply self$data was the matrix already or the df?
  # In R6 initialize: `self$terms = terms(response ~ ., data = self$data)` after `self$data` became matrix?
  # `terms.default` works on matrices.
  model_terms = tryCatch(terms(response ~ ., data = as.data.frame(data_mat)), error = function(e) NULL)
  
  # ------------------------------------------------------------------------
  # 2. Data Preparation
  # ------------------------------------------------------------------------
  
  X_concordant = NULL
  y_concordant = NULL
  treatment_concordant = NULL
  
  X_diffs_discordant = NULL
  y_diffs_discordant = NULL
  treatment_diffs_discordant = NULL
  
  # Internal variables for processing
  y = response
  X = data_mat
  X_model_matrix_col_names = colnames(X)
  
  if (length(y) <= ncol(X) + 5){
    stop("Not enough rows. Must be at least the number of covariates plus 5")
  }
  
  if (is.null(treatment)) {
    warning("Function did not include a treatment parameter. The order of the differences in the matched pairs will be arbitrary.")
  }
  
  matched_data = process_matched_pairs_cpp(
    strata = strata,
    y = y,
    X = X,
    treatment = treatment
  )
  
  X_concordant = matched_data$X_concordant
  y_concordant = matched_data$y_concordant
  treatment_concordant = matched_data$treatment_concordant
  strata_concordant = matched_data$strata_concordant
  X_diffs_discordant = matched_data$X_diffs_discordant
  y_diffs_discordant = matched_data$y_diffs_discordant
  treatment_diffs_discordant = matched_data$treatment_diffs_discordant
  
  # ------------------------------------------------------------------------
  # 3. Model Fitting
  # ------------------------------------------------------------------------
  
  discordant_model = NULL
  converged_discordant = NULL
  

  # --- Concordant Pairs / Reservoir Model ---
  
  concordant_model = NULL
  concordant_converged = NULL
  
  if (use_concordant_pairs) {
    # Check if we should warn about dropping concordant
    if (length(y_concordant) < ncol(X_concordant) + 5) {
      # This logic mimics R6 line 188. "droped_reservoir_concordant" there likely came from process_matched_pairs logic
      warning("There are not enough concordant pairs or reservoir entries. Those rows will be dropped and only the discordant pairs will be used.")
    } else {
      if (concordant_method == "GLM") {
        concordant_model = glm(y_concordant ~ treatment_concordant + X_concordant, family = "binomial")
        b_con = summary(concordant_model)$coefficients[,1]
        Sigma_con = pmin(vcov(concordant_model), 20)
        eps = 1e-6         #added for stabilibty
        Sigma_con = Sigma_con + diag(eps, nrow(Sigma_con))
        
        b_con = c(0, b_con[-c(1,2)])
        Sigma_con = Sigma_con[-1,-1]; Sigma_con[1,] = 0; Sigma_con[, 1] = 0; Sigma_con[1,1] = 20
      } else if (concordant_method == "GEE") {
        concordant_model = geeglm(
          y_concordant ~ treatment_concordant + X_concordant,
          id    = strata_concordant,
          family = binomial(link = "logit"),
          corstr = "exchangeable",
          data   = data.frame(y_concordant, treatment_concordant, X_concordant, strata_concordant)
        )
        
        b_con = summary(concordant_model)$coefficients[,1]
        Sigma_con = pmin(vcov(concordant_model), 20)
        eps = 1e-6         #added for stabilibty
        Sigma_con = Sigma_con + diag(eps, nrow(Sigma_con))
        
        b_con = c(0, b_con[-c(1,2)])
        Sigma_con = Sigma_con[-1,-1]; Sigma_con[1,] = 0; Sigma_con[, 1] = 0; Sigma_con[1,1] = 20
      } else if (concordant_method == "GLMM") {
        concordant_model = glmmTMB(
          y_concordant ~ treatment_concordant + X_concordant + (1 | strata_concordant),
          family = binomial(),
          data   = data.frame(y_concordant, treatment_concordant, X_concordant, strata_concordant)
        )
        
        b_con = summary(concordant_model)$coefficients$cond[,1]
        Sigma_con = pmin(vcov(concordant_model)$cond, 20)
        eps = 1e-6         #added for stabilibty
        Sigma_con = Sigma_con + diag(eps, nrow(Sigma_con))
        
        b_con = c(0, b_con[-c(1,2)])
        Sigma_con = Sigma_con[-1,-1]; Sigma_con[1,] = 0; Sigma_con[, 1] = 0; Sigma_con[1,1] = 20
      }
    }
  }

  # --- Discordant Pairs Model ---
  if (!droped_discordant) {
    # Fit fast model
    discordant_model = fast_conditional_logistic_regression_with_var_cpp(cbind(1, X_diffs_discordant), y_diffs_discordant, j = -1)
    
    converged_discordant = discordant_model$converged
    
    if (!discordant_model$converged) {
      # Fallback to GLM
      y_discord_1_0 = as.numeric(y_diffs_discordant == 1)
      coef_table_discordant = summary(suppressWarnings(glm(y_discord_1_0 ~ X_diffs_discordant, family = "binomial")))$coefficients
      
      if (coef_table_discordant[1, 2] < 1e4) {
        # GLM converged/stable
        discordant_model = list(b = coef_table_discordant[, 1],
                                ssq_b = coef_table_discordant[, 2]^2)
        converged_discordant = TRUE
        warning("Discordant pairs model is unstable but converged.")
      } else {
        # GLM Unstable
        if (try_to_stableize && use_concordant_pairs) {
          # Decompose discordant logic
          # Basically puts all original data into reservoir logic
          # R6 'decompose_discordant' basically reset reservoir to full data and cleared discordant
          # This seems to be because decomposing pairs -> just analyzing them as unmatched data
          
          # Logic from R6 decompose_discordant:
          # set X_reservoir = original data (which is data_mat/X)
          # set y_reservoir = original response (y)
          # set treatment_reservoir = original treatment
          
          X_concordant = X
          y_concordant = y
          treatment_concordant = treatment
          
          X_diffs_discordant = NULL
          y_diffs_discordant = NULL
          treatment_diffs_discordant = NULL
          
          droped_reservoir_concordant = TRUE # Wait, R6 set this to TRUE?
          # R6 line 177: self$droped_reservoir_concordant = TRUE.
          # The logic implies "we dropped the original reservoir/concordant plan and replaced it with everything?" 
          # Or "we dropped the idea of using separate concordant/reservoir and just used everything"?
          # Actually the flag name `droped_reservoir_concordant` is confusing in R6 context.
          
          discordant_model = NULL
          warning("Discordant pairs model did not converge. Decomposed the pairs and added them to the resevior.")
        } else {
          warning("Discordant pairs model did not converge.")
        }
      }
    }
  }
  
  # ------------------------------------------------------------------------
  # 4. Result Construction
  # ------------------------------------------------------------------------
  
  coefficients_and_variances = NULL
  
  if (is.null(concordant_model)) {
      if (!is.null(discordant_model)) {
          coefficients_and_variances = data.frame(
              coefs_discordant = discordant_model$b,
              std_err_discordant = sqrt(discordant_model$ssq_b)
          )
      }
  } else if (is.null(discordant_model)) {
       coefficients_and_variances = data.frame(
          coefs_reservoir_concordant = concordant_model$b,
          std_err_reservoir_concordant = sqrt(concordant_model$ssq_b)
       )
  } else {
      # Both converged
      w_star = concordant_model$ssq_b / (concordant_model$ssq_b + discordant_model$ssq_b)
      
      coefficients_and_variances = data.frame(
          coefs_discordant = discordant_model$b,
          std_err_discordant = sqrt(discordant_model$ssq_b),
          coefs_reservoir_concordant = concordant_model$b,
          std_err_reservoir_concordant = sqrt(concordant_model$ssq_b),
          coefs_combined = discordant_model$b * w_star + concordant_model$b * (1 - w_star),
          std_err_combined = sqrt(discordant_model$ssq_b * w_star)
      )
  }
  
  # Construct result object
  res = list(
      call = paste0("Call: clogitR(response = ", response_name, 
                     ", data = ", data_name, 
                     ", treatment = ", treatment_name, 
                     ", strata = ", strata_name, 
                     ", use_concordant_pairs = ", use_concordant_pairs,
                     ", try_to_stableize = ", try_to_stableize, ")"),
      coefficients_and_variances = coefficients_and_variances,
      converged_discordant = converged_discordant,
      concordant_converged = concordant_converged,
      use_concordant_pairs = use_concordant_pairs,
      num_discordant = length(y_diffs_discordant),
      num_concordant = length(y_concordant),
      X_model_matrix_col_names = X_model_matrix_col_names
  )
  
  class(res) = c("clogitR", "list")
  return(res)
}


#' @export
summary.clogitR = function(object, ...) {
  
  res = list()
  res$call = object$call
  res$converged_discordant = object$converged_discordant
  res$use_concordant_pairs = object$use_concordant_pairs
  res$concordant_converged = object$concordant_converged
  res$num_discordant = object$num_discordant
  res$num_concordant = object$num_concordant
  
  coeffs = object$coefficients_and_variances
  col_names = object$X_model_matrix_col_names
  
  # Note: The result vectors include Intercept (if present in X_diffs, which had 1 added in C++ calls)
  # In R6: `c("Intercept", private$X_model_matrix_col_names)`
  # fast_conditional... adds a column.
  # If data_mat has p columns, C++ gets p+1.
  # We should check dimensions.
  # Assuming standard behavior, the first coef is indeed intercept if added.
  # In R6 code: `fast_conditional_logistic_regression_with_var_cpp(cbind(1, ...))`
  # So yes, Intercept is index 1.
  
  row_names = c("Intercept", col_names)
  
  if (!is.null(coeffs$coefs_combined)) {
    mixed_coeffs = data.frame(
      coef = coeffs$coefs_combined,
      `exp(coef)` = exp(coeffs$coefs_combined),
      `Std. Error` = coeffs$std_err_combined,
      `z value` = coeffs$coefs_combined / coeffs$std_err_combined,
      pr =  2 * pnorm(-abs(coeffs$coefs_combined / coeffs$std_err_combined))
    )
    colnames(mixed_coeffs)[4] = "Pr(>|z|)"
    rownames(mixed_coeffs) = row_names
    res$mixed_coefficients = mixed_coeffs
  }
  
  if (!is.null(coeffs$coefs_discordant)) {
    discordant_coeffs = data.frame(
      coef = coeffs$coefs_discordant,
      `exp(coef)` = exp(coeffs$coefs_discordant),
      `Std. Error` = coeffs$std_err_discordant,
      `z value` = coeffs$coefs_discordant / coeffs$std_err_discordant,
      pr =  2 * pnorm(-abs(coeffs$coefs_discordant / coeffs$std_err_discordant))
    )
    colnames(discordant_coeffs)[4] = "Pr(>|z|)"
    rownames(discordant_coeffs) = row_names
    res$discordant_coefficients = discordant_coeffs
  }
  
  if (!is.null(coeffs$coefs_reservoir_concordant) && object$use_concordant_pairs) {
    concordant_coeffs = data.frame(
      coef = coeffs$coefs_reservoir_concordant,
      `exp(coef)` = exp(coeffs$coefs_reservoir_concordant),
      `Std. Error` = coeffs$std_err_reservoir_concordant,
      `z value` = coeffs$coefs_reservoir_concordant / coeffs$std_err_reservoir_concordant,
      pr =  2 * pnorm(-abs(coeffs$coefs_reservoir_concordant / coeffs$std_err_reservoir_concordant))
    )
    colnames(concordant_coeffs)[4] = "Pr(>|z|)"
    rownames(concordant_coeffs) = row_names
    res$concordant_coefficients = concordant_coeffs
  }
  
  class(res) = "summary.clogitR"
  res
}

#' @export
print.summary.clogitR = function(x, ...) {
  
  cat(x$call, "\n\n")
  
  cat("Discordant model converged =", x$converged_discordant)
  if (x$use_concordant_pairs) { cat(" | Concordant and reservoir model converged =", x$concordant_converged, '\n\n') } else {cat("\n\n")}
  
  if (!is.null(x$mixed_coefficients)) {
    cat("Mixed Coefficients:\n")
    printCoefmat(x$mixed_coefficients, digits = 4, signif.stars = TRUE)
    cat("\n")
  }
  
  if (!is.null(x$discordant_coefficients)) {
    cat("Discordant Coefficients:\n")
    printCoefmat(x$discordant_coefficients, digits = 4, signif.stars = TRUE)
    cat("\n")
  }
  
  if (!is.null(x$concordant_coefficients)) {
    cat("Concordant and Reservoir Coefficients:\n")
    printCoefmat(x$concordant_coefficients, digits = 4, signif.stars = TRUE)
    cat("\n")
  }
  
  cat("Number of discordant pairs =", x$num_discordant,
      ", Number of condordant and reservoir entries =", x$num_concordant)
  
  invisible(x)
}

#' @export
process_matched_pairs = function(
    strata = NULL,
    y = NULL,
    X = NULL,
    treatment = NULL
) {
  process_matched_pairs_cpp(
    strata = strata,
    y = y,
    X = X,
    treatment = treatment
  )
}