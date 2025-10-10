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
    
    #' @field formula The model formula supplied by the user.
    formula = NULL,
    
    #' @field data_name The name of the dataset provided to the model.
    data_name = NULL,
    
    #' @field data The dataset provided to the model.
    data = NULL,
    
    #' @field model_frame The model frame generated from the formula and data.
    model_frame = NULL,
    
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
    
    #' @field X_diffs_discordant Difference of X rows in the strata. If treatment was indicated in the formula the difference is T - C.
    X_diffs_discordant = NULL,
    
    #' @field y_diffs_discordant Difference of y values in the strata. If treatment was indicated in the formula the difference is T - C.
    y_diffs_discordant = NULL,
    
    #' @field coefficients_and_variances Table containing the coefficients and variances of each covariate.
    coefficients_and_variances = NULL,
    
    #' @field droped_reservoir_concordant Flag indicating whether the reservoir and concordant data was doped.
    droped_reservoir_concordant = NULL,
    
    #' @field droped_discordant Flag indicating whether the discordant data was doped.
    droped_discordant = NULL,
    
    #' @field converged_reservoir_concordant Flag indicating whether the reservoir and concordant model converged.
    converged_reservoir_concordant = NULL,
    
    #' @field converged_discordant Flag indicating whether the discordant model converged.
    converged_discordant = NULL,
    
    #' Initialize a new clogitRClass object
    #'
    #' @param formula A formula specifying the model.
    #' @param data A data.frame or data.table containing the variables.
    #' @param use_concordant_pairs Flag indicating whether to include the concordant pairs and reservoir in the model.
    #' @param try_to_stableize Flag indicating whether the model should decompose the discordant pairs to increase stability if the model is unstable.
    #' @return A new `clogitRClass` object.
    initialize = function(formula,
                          data,
                          use_concordant_pairs = TRUE,
                          try_to_stableize = TRUE
                          ) {
      assertFormula(formula)
      
      if (!is.null(data)) {
        assertDataFrame(data)
        self$formula = formula
        self$data = data
        self$model_frame = model.frame(formula, data)
        self$terms = terms(formula, data)
        
        # Extract response
        self$y = model.response(self$model_frame)
        
        # Extract predictors (matrix)
        self$X = model.matrix(delete.response(self$terms), self$model_frame)
        
      } 
      
      
      
      assertDataFrame(data)
      assertFlag(use_concordant_pairs)
      assertFlag(try_to_stableize)
      
      self$formula = formula
      self$data_name <- deparse(substitute(data))
      self$data = data
      self$use_concordant_pairs = use_concordant_pairs
      self$try_to_stableize = try_to_stableize
      
      # Build model frame and terms
      self$model_frame = model.frame(formula, data)
      self$terms = terms(formula, data = data)
      
      slef$fit()
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
        discordant_model = fast_conditional_logistic_regression_with_var_cpp(self$X_diffs, self$y_diffs, j = -1)

        self$converged_discordant = discordant_model$converged
        if (!discordant_model$converged){
          #If the fast model did not converge, attempt to a standard glm model (slower)
          y_discord_1_0 = as.numeric(self$y_diffs == 1)
          coef_table_discordant = summary(suppressWarnings(glm(y_discord_1_0 ~ 0 + self$X_diffs, family = "binomial")))$coefficients
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
              self$droped_reservoir_concordant = FALSE
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
        y_1_n1_reservoir = 2 * self$y_reservoir - 1
        concordant_model = fast_conditional_logistic_regression_with_var_cpp(cbind(1, self$X_reservoir), self$y_1_n1_reservoir, j = -1)
        self$converged_reservoir_concordant = concordant_model$converged
        if (!concordant_model$converged){
          #If the fast model did not converge, attempt to a standard glm model (slower)
          coef_table_concordant = summary(suppressWarnings(glm(y_reservoir ~ self$X_diffs, family = "binomial")))$coefficients
          
          if (coef_table_concordant[2,2] < 1e4){
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
        #the concordant model has an intercept so we have to get rid of it
        concordant_model$b = concordant_model$b[2:length(concordant_model$b)]
        concordant_model$ssq_b = concordant_model$ssq_b[2:length(concordant_model$ssq_b)]
      }
      
      w_star = concordant_model$ssq_b / (concordant_model$ssq_b + discordant_model$ssq_b)
      
      self$coefficients_and_variances = data.frame(
        coefs_discordant = discordant_model$b,
        std_err_discordant = sqrt(discordant_model$ssq_b),
        coefs_reservoir_concordant = concordant_model$b,
        std_err_reservoir_concordant = sqrt(concordant_model$ssq_b),
        coefs_combined = discordant_model$b * w_star + concordant_model$b * (1 - w_star),
        std_err_reservoir_combined = sqrt(discordant_model$ssq_b * w_star)
        )
    },
    
    #' Summarize the model
    #'
    #' @description
    #' Placeholder summary function.
    summary = function() {
      print("Call:")
      call = cat("clogitR(", self$formula, ", data = ", self$data_name, ", use_concordant_pairs = ", self$use_concordant_pairs, ", try_to_stableize = ", self$try_to_stableize, ")")
      cat(call, '\n')
      
      if (self$converged_discordant) { cat("Discordant model converged = ", self$converged_discordant, '\n') }
      if (self$use_concordant_pairs & self$converged_reservoir_concordant) { cat("Concordant and reservoir model converged = ", self$converged_reservoir_concordant, '\n') }
      
      
      if (self$converged_discordant & self$use_concordant_pairs & self$converged_reservoir_concordant) {
        
        mixed_coeffs = data.frame(
          Estimate = self$coefficients_and_variances$coefs_combined,
          `Std. Error` = self$coefficients_and_variances$std_err_reservoir_combined,
          `z value` = self$coefficients_and_variances$coefs_combined / self$coefficients_and_variances$std_err_reservoir_combined,
          `Pr(>|z|)` =  2 * pnorm(-abs(self$coefficients_and_variances$coefs_combined / self$coefficients_and_variances$std_err_reservoir_combined))
        )
        rownames(concordant_coeffs) = private$X_model_matrix_col_names
        print("Concordant and Reservoir Coefficients:")
        printCoefmat(concordant_coeffs, digits = 4, signif.stars = TRUE)
        
      } else if (self$converged_discordant) {
        
        discordant_coeffs = data.frame(
          Estimate = self$coefficients_and_variances$coefs_discordant,
          `Std. Error` = self$coefficients_and_variances$std_err_discordant,
          `z value` = self$coefficients_and_variances$coefs_discordant / self$coefficients_and_variances$std_err_discordant,
          `Pr(>|z|)` =  2 * pnorm(-abs(self$coefficients_and_variances$coefs_discordant / self$coefficients_and_variances$std_err_discordant))
        )
        rownames(discordant_coeffs) = private$X_model_matrix_col_names
        print("Discordant Coefficients:")
        printCoefmat(discordant_coeffs, digits = 4, signif.stars = TRUE)
        
      } else {
        
        concordant_coeffs = data.frame(
          Estimate = self$coefficients_and_variances$coefs_reservoir_concordant,
          `Std. Error` = self$coefficients_and_variances$std_err_reservoir_concordant,
          `z value` = self$coefficients_and_variances$coefs_reservoir_concordant / self$coefficients_and_variances$std_err_reservoir_concordant,
          `Pr(>|z|)` =  2 * pnorm(-abs(self$coefficients_and_variances$coefs_reservoir_concordant / self$coefficients_and_variances$std_err_reservoir_concordant))
        )
        rownames(concordant_coeffs) = private$X_model_matrix_col_names
        print("Concordant and Reservoir Coefficients:")
        printCoefmat(concordant_coeffs, digits = 4, signif.stars = TRUE)
      }
      
      
      #cat("\n Bryar score: ", ) #TODO
      
      list(call = call,
           mixed_coefficients = mixed_coeffs,
           discordant_coefficients = discordant_coeffs,
           concordant_coefficients = concordant_coeffs
           )
    }
  ),
  
  private = list(
    X_model_matrix_col_names = list(),
    
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
      
      self$X_reservoir = X
      self$y_reservoir = y
    },
    
    prepare_data = function() {
      # Extract response and predictors
      y = model.response(self$model_frame)
      X = model.matrix(update(self$terms, . ~ . -1), self$model_frame)

      if (length(y) <= ncol(X) + 5){
        stop("Not enough rows. Must be at least the number of covariates plus 5")
      }
      
      # Find treatment
      treatment_var = attr(self$terms, "specials")$treatment
      if (!is.null(treatment_var)) {
        treatment_name = all.vars(attr(self$terms, "variables")[[treatment_var + 1]])
        treatment = self$model_frame[[treatment_name]]
        if (!all(treatment %in% c(0,1))) {
          stop("Treatment must be binary 0 or 1.")
        }
        if (!(treatment_name %in% colnames(X))) {
          stop("Treatment variable not found in model matrix.")
        }
        X = X[, c(treatment_name, setdiff(colnames(X), treatment_name)), drop = FALSE]
      } else {
        warning("Formula did not include a treatment() term. The order of the differences in the matched pairs will be arbitrary.")
      }
      private$X_model_matrix_col_names = colnames(X)
      
      # Find strata
      strata_var = attr(self$terms, "specials")$strata
      if (is.null(strata_var)) {
        self$X_reservoir = X
        self$y_reservoir = y
        self$droped_discordant = TRUE
        warning("Formula did not include a strata() term. All the data was put in the reservoir.")
        return(NULL)
      }
      strata_name = all.vars(attr(self$terms, "variables")[[strata_var + 1]])
      strata = self$model_frame[[strata_name]]
      
      #remove the strata from X
      if (strata_name %in% colnames(X)) {
        X = X[, setdiff(colnames(X), strata_name), drop = FALSE]
      }
      
      # Initialize storage
      reservoir_X = list()
      reservoir_y = list()
      diffs_X = list()
      diffs_y = list()
      
      # Split reservoir vs matched
      reservoir_idx = which(strata == 0)
      matched_idx = which(strata != 0)
      
      #checking for valid sizes
      #check for enough discordant pairs
      if((length(matched_idx) / 2) <= p + 5){
        self$X_reservoir = X
        self$y_reservoir = y
        self$droped_discordant = TRUE
        warning("There are not enough discordant pairs. All the data was put in the reservoir.")
        return(NULL)
      }
      
      #check for inclusion of the reservoir
      if(!self$use_concordant_pairs){
        reservoir_idx = NULL
      }
      
      #check for enough concordant and reservior
      if(length(reservoir_idx) <= p + 5){
        reservoir_idx = NULL
        self$droped_reservoir_concordant = TRUE
        #Throw the warning for this only when we are sure we are not adding the discordant pairs to the reservoir
        #warning("There are not enough concordant pairs or reservoir entries. Those rows will be dropped and only the discordant pairs will be used.")
      }
      
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
      self$X_reservoir_concordant = if (length(reservoir_X) > 0) do.call(rbind, reservoir_X) else NULL
      self$y_reservoir_concordant = if (length(reservoir_y) > 0) unlist(reservoir_y) else NULL
      self$X_diffs_discordant = if (length(diffs_X) > 0) do.call(rbind, diffs_X) else NULL
      self$y_diffs_discordant = if (length(diffs_y) > 0) unlist(diffs_y) else NULL
    }
  )
)