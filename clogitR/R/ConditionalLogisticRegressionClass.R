#' Conditional Logistic Regression
#'
#' @description
#' `clogitRClass` is an R6 class for fitting conditional logistic regression
#' models using stratified case-control data. This is a skeleton implementation
#' containing structure and documentation, but no estimation routine.
#'
#' @details
#' The class accepts a formula interface (`y ~ x1 + x2 + strata(stratum)`)
#' and a `data.frame` or `data.table`. Internally, it parses the formula,
#' extracts variables, and prepares the data for model fitting.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   y = c(1, 0, 1, 0),
#'   x = c(2.3, 1.5, 3.1, 2.9),
#'   stratum = c(1, 1, 2, 2)
#' )
#'
#' fit <- clogitR$new(y ~ x + strata(stratum), data = data)
#' fit$summary()
#' }
#'
clogitRClass = R6::R6Class(
  classname = "clogitRClass",
  
  public = list(
    
    #' @field response Binary vector indicating the response of each subject.
    response = NULL,
    
    #' @field strata Numeric vector indexing the matched pairs and the reservoir.
    strata = NULL,
    
    #' @field treatment Binary vector indicating weather a subject is test or control.
    treatment = NULL,
    
    #' @field data The dataset provided to the model.
    data = NULL,
    
    #' @field terms Model terms object from the formula.
    terms = NULL,
    
    #' @field use_concordant_pairs Indicator to fit the concordant pairs.
    use_concordant_pairs = NULL,
    
    #' @field try_to_stableize Indicator to decompose the discordant pairs the pairs if the model is unstable.
    try_to_stableize = NULL,
    
    #' @field X_reservoir_concordant X rows marked with 0 in the strata or where part of the concordat pairs.
    X_reservoir_concordant = NULL,
    
    #' @field y_reservoir_concordant y values marked with 0 in the strata or where part of the concordat pairs.
    y_reservoir_concordant = NULL,
    
    #' @field treatment_reservoir_concordant Treatment values marked with 0 in the strata or where part of the concordat pairs.
    treatment_reservoir_concordant = NULL,
    
    #' @field X_diffs_discordant Difference of X rows in the strata. If treatment was indicated in the formula the difference is T - C.
    X_diffs_discordant = NULL,
    
    #' @field y_diffs_discordant Difference of y values in the strata. If treatment was indicated in the formula the difference is T - C.
    y_diffs_discordant = NULL,
    
    #' @field treatment_diffs_discordant Difference of treatment values in the strata. This should only ever be `NULL` or all 1's.
    treatment_diffs_discordant = NULL,
    
    #' @field coefficients_and_variances Table containing the coefficients and variances of each covariate.
    coefficients_and_variances = NULL,
    
    #' @field droped_reservoir_concordant Flag indicating whether the reservoir and concordant data was doped.
    droped_reservoir_concordant = FALSE,
    
    #' @field droped_discordant Flag indicating whether the discordant data was doped.
    droped_discordant = FALSE,
    
    #' @field converged_reservoir_concordant Flag indicating whether the reservoir and concordant model converged.
    converged_reservoir_concordant = NULL,
    
    #' @field converged_discordant Flag indicating whether the discordant model converged.
    converged_discordant = NULL,
    
    #' Initialize a new clogitRClass object
    #'
    #' @param response A binary (0,1) vector containing the response of each subject.
    #' @param data A data.frame, data.table, or model.matrix containing the variables.
    #' @param treatment A binary (0,1) vector containing weather each subject is treatment or control.
    #' @param strata Numeric vector indexing the matched pairs and the reservoir. Entrees indexed with zero, are added to the reservoir. Entrees indexed with the same number are part of the same strata.
    #' @param use_concordant_pairs Flag indicating whether to include the concordant pairs and reservoir in the model.
    #' @param try_to_stableize Flag indicating whether the model should decompose the discordant pairs to increase stability if the model is unstable.
    #' @return A new `clogitRClass` object.
    initialize = function(response = NULL,
                          data = NULL,
                          treatment = NULL,
                          strata = NULL,
                          use_concordant_pairs = TRUE,
                          try_to_stableize = TRUE
                          ) {
      if (test_data_frame(data, types = c("numeric", "integer", "factor", "logical"))) {
        self$data = model.matrix(~0+., data = data)
        
        # Case 2: already a model matrix (numeric matrix)
      } else if (test_matrix(data, mode = "numeric")) {
        if (sd(data[, 1]) == 0) { data[, 1] = NULL }
        self$data = as.matrix(data)
        
        # Otherwise: invalid type
      } else {
        assert(
          check_data_frame(data),
          check_matrix(data, mode = "numeric")
        )
      }
      n = nrow(data)
      private$data_name = deparse(substitute(data))
      
      assertFlag(use_concordant_pairs)
      self$use_concordant_pairs = use_concordant_pairs
      
      assertFlag(try_to_stableize)
      self$try_to_stableize = try_to_stableize
      
      assertNumeric(response, lower = 0, upper = 1, any.missing = FALSE, len = n)
      self$response = response
      private$response_name = deparse(substitute(response))
      
      if (!is.null(treatment)) {
        assertNumeric(treatment, lower = 0, upper = 1, any.missing = FALSE, len = n)
        if (!all(treatment %in% c(0,1))) { stop("Treatment must be binary 0 or 1.") }
        self$treatment = treatment
        private$treatment_name = deparse(substitute(treatment))
      }
      
      if (!is.null(strata)) {
        assertNumeric(strata, any.missing = FALSE, len = n)
        self$strata = strata
        private$strata_name  = deparse(substitute(strata))
      }
      
      
      # Build model frame and terms
      self$terms = terms(response ~ ., data = self$data)
      
      self$fit()
    },
    
    #' Fit the model
    #'
    #' @description
    #' Placeholder function for estimation. Currently not implemented.
    fit = function() {
      
      if (is.null(self$X_reservoir_concordant) && is.null(self$X_diffs_discordant)) {
        private$prepare_data()
      }
      
      if (!self$droped_discordant){
        #Fit fast model
        discordant_model = fast_conditional_logistic_regression_with_var_cpp(cbind(1, self$X_diffs_discordant), self$y_diffs_discordant, j = -1)

        self$converged_discordant = discordant_model$converged
        if (!discordant_model$converged){
          #If the fast model did not converge, attempt to a standard glm model (slower)
          y_discord_1_0 = as.numeric(self$y_diffs == 1)
          coef_table_discordant = summary(suppressWarnings(glm(y_discord_1_0 ~ self$X_diffs_discordant, family = "binomial")))$coefficients
          if (coef_table_discordant[1,2] < 1e4){
            #If the glm model is good, just use that
            discordant_model = list(b = coef_table_discordant[,1],
                                    ssq_b = coef_table_discordant[,2]^2)
            self$converged_discordant = TRUE
            warning("Discordant pairs model is unstable but converged.")
          } else {
            #If the glm model is also bad we give up on the discordant pairs (there just weren't enough)
            if (self$try_to_stableize && self$use_concordant_pairs) {
              private$decompose_discordant()
              self$droped_reservoir_concordant = TRUE
              discordant_model = NULL
              warning("Discordant pairs model did not converge. Decomposed the pairs and added them to the resevior.")
            } else {
              warning("Discordant pairs model did not converge.")
            }
          }
        }
      }
      
      #only throw the warning of dropping the concordant pairs if we are sure we are not adding the discordant pairs to it
      if (self$droped_reservoir_concordant && self$use_concordant_pairs && self$converged_discordant) {
        warning("There are not enough concordant pairs or reservoir entries. Those rows will be dropped and only the discordant pairs will be used.")
      }
      
      if (self$use_concordant_pairs && !self$droped_reservoir_concordant) {
        #first try the fast model
        y_1_n1_reservoir = 2 * self$y_reservoir_concordant - 1
        
        concordant_model =
          if (!is.null(self$treatment_reservoir_concordant)) {
            fast_conditional_logistic_regression_with_var_cpp(cbind(1, self$treatment_reservoir_concordant, self$X_reservoir_concordant), y_1_n1_reservoir, j = -1)
          } else {
            fast_conditional_logistic_regression_with_var_cpp(cbind(1, self$X_reservoir_concordant), y_1_n1_reservoir, j = -1)
          }
        
        self$converged_reservoir_concordant = concordant_model$converged
        if (!concordant_model$converged){
          #If the fast model did not converge, attempt to a standard glm model (slower)
          coef_table_concordant =
            if (!is.null(self$treatment_reservoir_concordant)) {
              summary(suppressWarnings(glm(self$y_reservoir_concordant ~ self$treatment_reservoir_concordant + self$X_reservoir_concordant, family = "binomial")))$coefficients
            } else {
              summary(suppressWarnings(glm(self$y_reservoir_concordant ~ self$X_reservoir_concordant, family = "binomial")))$coefficients
            }
          
          if ((!is.null(self$treatment_reservoir_concordant) & coef_table_concordant[2,2] < 1e4) || (is.null(self$treatment_reservoir_concordant) & coef_table_concordant[1,2] < 1e4)){
            #If the glm model is good, just use that
            discordant_model = list(b = coef_table_concordant[,1],
                                    ssq_b = coef_table_concordant[,2]^2)
            self$converged_reservoir_concordant = TRUE
            warning("Concordant pairs model is unstable but converged.")
          } else {
            #If the glm model is also bad we give up on the concordant pairs (there just weren't enough)
            concordant_model = NULL
            warning("Discordant pairs model did not converge. Only the discordant pairs model will be used.")
          }
        }
        #if there is a treatment, then the concordant mod has an intersept and a w col. we have to get rid of it to make everything align
        if (!is.null(self$treatment_reservoir_concordant)) {
          concordant_model$b = concordant_model$b[2:length(concordant_model$b)]
          concordant_model$ssq_b = concordant_model$ssq_b[2:length(concordant_model$ssq_b)]
        }
      }
      
      if (is.null(concordant_model)){
        self$coefficients_and_variances = data.frame(
          coefs_discordant = discordant_model$b,
          std_err_discordant = sqrt(discordant_model$ssq_b)
        )
      } else if (is.null(discordant_model)) {
        self$coefficients_and_variances = data.frame(
          coefs_reservoir_concordant = concordant_model$b,
          std_err_reservoir_concordant = sqrt(concordant_model$ssq_b)
        )
      } else {
        w_star = concordant_model$ssq_b / (concordant_model$ssq_b + discordant_model$ssq_b)
        
        self$coefficients_and_variances = data.frame(
          coefs_discordant = discordant_model$b,
          std_err_discordant = sqrt(discordant_model$ssq_b),
          coefs_reservoir_concordant = concordant_model$b,
          std_err_reservoir_concordant = sqrt(concordant_model$ssq_b),
          coefs_combined = discordant_model$b * w_star + concordant_model$b * (1 - w_star),
          std_err_combined = sqrt(discordant_model$ssq_b * w_star)
        )
      }
    },
    
    #' Summarize the model
    #'
    #' @description
    #' Placeholder summary function.
    summary_data = function() {
      ret = list()
      ret$call = paste0("Call: clogitR(response = ", private$response_name, 
                     ", data = ", private$data_name, 
                     ", treatment = ", private$treatment_name, 
                     ", strata = ", private$strata_name, 
                     ", use_concordant_pairs = ", self$use_concordant_pairs,
                     ", try_to_stableize = ", self$try_to_stableize, ")")
      ret$converged_discordant = self$converged_discordant
      ret$use_concordant_pairs = self$use_concordant_pairs
      ret$converged_reservoir_concordant = self$converged_reservoir_concordant
      
      
      if (!is.null(self$coefficients_and_variances$coefs_combined)) {
        
        mixed_coeffs = data.frame(
          coef = self$coefficients_and_variances$coefs_combined,
          `exp(coef)` = exp(self$coefficients_and_variances$coefs_combined),
          `Std. Error` = self$coefficients_and_variances$std_err_combined,
          `z value` = self$coefficients_and_variances$coefs_combined / self$coefficients_and_variances$std_err_combined,
          pr =  2 * pnorm(-abs(self$coefficients_and_variances$coefs_combined / self$coefficients_and_variances$std_err_combined))
        )
        colnames(mixed_coeffs)[4] = "Pr(>|z|)"
        rownames(mixed_coeffs) = c("Intercept", private$X_model_matrix_col_names)
        ret$mixed_coefficients = mixed_coeffs
      }
      if (!is.null(self$coefficients_and_variances$coefs_discordant)) {
        
        discordant_coeffs = data.frame(
          coef = self$coefficients_and_variances$coefs_discordant,
          `exp(coef)` = exp(self$coefficients_and_variances$coefs_discordant),
          `Std. Error` = self$coefficients_and_variances$std_err_discordant,
          `z value` = self$coefficients_and_variances$coefs_discordant / self$coefficients_and_variances$std_err_discordant,
          pr =  2 * pnorm(-abs(self$coefficients_and_variances$coefs_discordant / self$coefficients_and_variances$std_err_discordant))
        )
        colnames(discordant_coeffs)[4] = "Pr(>|z|)"
        rownames(discordant_coeffs) = c("Intercept", private$X_model_matrix_col_names)
        ret$discordant_coefficients = discordant_coeffs
      } 
      if (!is.null(self$coefficients_and_variances$coefs_reservoir_concordant) & self$use_concordant_pairs) {
        
        concordant_coeffs = data.frame(
          coef = self$coefficients_and_variances$coefs_reservoir_concordant,
          `exp(coef)` = exp(self$coefficients_and_variances$coefs_reservoir_concordant),
          `Std. Error` = self$coefficients_and_variances$std_err_reservoir_concordant,
          `z value` = self$coefficients_and_variances$coefs_reservoir_concordant / self$coefficients_and_variances$std_err_reservoir_concordant,
          pr =  2 * pnorm(-abs(self$coefficients_and_variances$coefs_reservoir_concordant / self$coefficients_and_variances$std_err_reservoir_concordant))
        )
        colnames(concordant_coeffs)[4] = "Pr(>|z|)"
        rownames(concordant_coeffs) = c("Intercept", private$X_model_matrix_col_names)
        ret$concordant_coefficients = concordant_coeffs
      }
      
      ret$num_discordant = length(self$y_diffs_discordant)
      ret$num_concordant = length(self$y_reservoir_concordant)
      
      class(ret) = "summary.clogitR"
      ret
    }
  ),
  
  private = list(
    X_model_matrix_col_names = list(),
    data_name = NULL,
    response_name = NULL,
    treatment_name = NULL,
    strata_name = NULL,
    
    decompose_discordant = function() {
      y = model.response(self$model_frame)
      X = model.matrix(update(self$terms, . ~ . -1), self$model_frame)
      
      # Find treatment
      treatment_var = attr(self$terms, "specials")$treatment
      if (!is.null(treatment_var)) {
        treatment_name = all.vars(attr(self$terms, "variables")[[treatment_var + 1]])
        treatment = self$model_frame[[treatment_name]]
        X = X[, c(treatment_name, setdiff(colnames(X), treatment_name)), drop = FALSE]
      }
      print("decompose")
      self$X_reservoir_concordant = self$data
      self$y_reservoir_concordant = self$response
      if (!is.null(treatment)) {
        self$treatment_reservoir_concordant = treatment
      }
      self$X_diffs_discordant = NULL
      self$y_diffs_discordant = NULL
      self$treatment_diffs_discordant = NULL
    },
    
    prepare_data = function() {
      # Extract response and predictors
      y = self$response
      X = self$data

      if (length(y) <= ncol(X) + 5){
        stop("Not enough rows. Must be at least the number of covariates plus 5")
      }
      
      # Find treatment
      treatment = self$treatment
      if (is.null(treatment)) {
        warning("Function did not include a treatment parameter. The order of the differences in the matched pairs will be arbitrary.")
      }
      private$X_model_matrix_col_names = colnames(X)
      
      # Find strata
      strata = self$strata
      if (is.null(strata)) {
        print("added")
        self$X_reservoir_concordant = X
        self$y_reservoir_concordant = y
        if (!is.null(treatment)) {
          self$treatment_reservoir_concordant = treatment
        }
        self$droped_discordant = TRUE
        warning("Function did not include a strata parameter. All the data was put in the reservoir.")
        return(NULL)
      }
      
      matched_data =
        process_matched_pairs_cpp(
        strata = strata,
        y = y,
        X = X,
        treatment = treatment
      )
      
      self$X_reservoir_concordant =         matched_data$X_reservoir_concordant
      self$y_reservoir_concordant =         matched_data$y_reservoir_concordant
      self$treatment_reservoir_concordant = matched_data$treatment_reservoir_concordant
      self$X_diffs_discordant =             matched_data$X_diffs_discordant
      self$y_diffs_discordant =             matched_data$y_diffs_discordant
      self$treatment_diffs_discordant =     matched_data$treatment_diffs_discordant
      self$droped_discordant =              matched_data$dropped_discordant
      self$droped_reservoir_concordant =    matched_data$dropped_reservoir_concordant
    }
  )
)