#' Conditional Logistic Regression
#'
#' @description
#' `clogitR` is an R6 class for fitting conditional logistic regression
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
clogitR = R6::R6Class(
  classname = "clogitR",
  
  public = list(
    
    #' @field formula The model formula supplied by the user.
    formula = NULL,
    
    #' @field data The dataset provided to the model.
    data = NULL,
    
    #' @field model_frame The model frame generated from the formula and data.
    model_frame = NULL,
    
    #' @field terms Model terms object from the formula.
    terms = NULL,
    
    #' @field use_concordant_pairs Indicator to fit the concordant pairs.
    use_concordant_pairs = NULL,
    
    #' @field X_reservoir X rows marked with 0 in the strata or where part of the concordat pairs.
    X_reservoir = NULL,
    
    #' @field y_reservoir y values marked with 0 in the strata or where part of the concordat pairs.
    y_reservoir = NULL,
    
    #' @field X_diffs Difference of X rows in the strata. If treatment was indicated in the formula the difference is T - C.
    X_diffs = NULL,
    
    #' @field y_diffs Difference of y values in the strata. If treatment was indicated in the formula the difference is T - C.
    y_diffs = NULL,
    
    #' Initialize a new clogitR object
    #'
    #' @param formula A formula specifying the model.
    #' @param data A data.frame or data.table containing the variables.
    #' @return A new `clogitR` object.
    initialize = function(formula,
                          data,
                          use_concordant_pairs = TRUE
                          ) {
      assertFormula(formula)
      assertDataFrame(data)
      assertFlag(use_concordant_pairs)
      
      self$formula = formula
      self$data = data
      self$use_concordant_pairs = use_concordant_pairs
      
      # Build model frame and terms
      self$model_frame = model.frame(formula, data)
      self$terms = terms(formula, data = data)
    },
    
    #' Fit the model
    #'
    #' @description
    #' Placeholder function for estimation. Currently not implemented.
    fit = function() {
      discordant_model = fast_conditional_logistic_regression_with_var_cpp(self$X_diffs,
                                                                           self$y_diffs,
                                                                           j = -1)
      if (!self$use_concordant_pairs) {
        private$beta_hat = discordant_model$b
        private$sse_beta_hat = discordant_model$ssq_b
      } else {
        y_1_n1_reservoir = 2 * y_reservoir - 1
        concordant_model = fast_conditional_logistic_regression_with_var_cpp(self$X_reservoir,
                                                                             self$y_1_n1_reservoir,
                                                                             j = -1)
        for (i in 1:length(discordant_model$b)) {
          
        }
      }
    },
    
    #' Summarize the model
    #'
    #' @description
    #' Placeholder summary function.
    summary = function() {
      cat("Conditional Logistic Regression Model\n")
      cat("Formula: ", deparse(self$formula), "\n")
      cat("Number of observations: ", nrow(self$model_frame), "\n")
      cat("Terms: ", paste(attr(self$terms, "term.labels"), collapse = ", "), "\n")
    }
  ),
  
  private = list(
    beta_hat = list(),
    sse_beta_hat = list(),
    
    prepare_data = function() {
      # Extract response and predictors
      y = model.response(self$model_frame)
      X = model.matrix(self$terms, self$model_frame)
      
      # Find strata
      strata_var = attr(self$terms, "specials")$strata
      if (is.null(strata_var)) {
        stop("Formula must include a strata() term.")
      }
      strata_name = all.vars(attr(self$terms, "variables")[[strata_var + 1]])
      strata = self$model_frame[[strata_name]]
      
      # Find treatment
      treatment_var = attr(self$terms, "specials")$treatment
      if (!is.null(treatment_var)) {
        treatment_name = all.vars(attr(self$terms, "variables")[[treatment_var + 1]])
        treatment = self$model_frame[[treatment_name]]
      }
      
      
      # Initialize storage
      reservoir_X = list()
      reservoir_y = list()
      diffs_X = list()
      diffs_y = list()
      
      # Split reservoir vs matched
      reservoir_idx = which(strata == 0)
      matched_idx = which(strata != 0)
      
      if (length(reservoir_idx) > 0) {
        reservoir_X[[length(reservoir_X) + 1]] = X[reservoir_idx, , drop = FALSE]
        reservoir_y[[length(reservoir_y) + 1]] = y[reservoir_idx]
      }
      
      # Handle matched pairs
      matched_list = split(matched_idx, strata[matched_idx])
      
      diffs_X = list()
      diffs_y = list()
      
      for (pair in matched_list) {
        if (length(pair) != 2) {
          stop("Each nonzero stratum must have exactly 2 rows.")
        }
        i = pair[1]
        j = pair[2]
        
        yi = y[i]
        yj = y[j]
        
        if(!is.null(treatment)){
          wi = treatment[i]
          wj = treatment[j]
        }
        
        if (yi == yj) {
          # Concordant pair → move both to reservoir
          reservoir_X[[length(reservoir_X) + 1]] = X[c(i, j), , drop = FALSE]
          reservoir_y[[length(reservoir_y) + 1]] = y[c(i, j)]
        } else {
          # Discordant → compute diff (case - control)
          if (!is.null(treatment)){
            if (wi == 1){
              diffs_X[[length(diffs_X) + 1]] = X[i, ] - X[j, ]
              diffs_y[[length(diffs_y) + 1]] = y[i] - y[j]
            } else {
              diffs_X[[length(diffs_X) + 1]] = X[j, ] - X[i, ]
              diffs_y[[length(diffs_y) + 1]] = y[j] - y[i]
            }
          } else {
            diffs_X[[length(diffs_X) + 1]] = X[i, ] - X[j, ]
            diffs_y[[length(diffs_y) + 1]] = y[i] - y[j]
          }
        }
      }
      
      # Store results
      self$X_reservoir = if (length(reservoir_X) > 0) do.call(rbind, reservoir_X) else NULL
      self$y_reservoir = if (length(reservoir_y) > 0) unlist(reservoir_y) else NULL
      self$X_diffs = if (length(diffs_X) > 0) do.call(rbind, diffs_X) else NULL
      self$y_diffs = if (length(diffs_y) > 0) unlist(diffs_y) else NULL
    }
  )
)


