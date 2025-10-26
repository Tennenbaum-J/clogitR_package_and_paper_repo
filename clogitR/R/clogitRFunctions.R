#' Initialize a new clogitR model
#'
#' @param response A binary (0,1) vector containing the response of each subject. HI
#' @param data A data.frame, data.table, or model.matrix containing the variables.
#' @param treatment A binary (0,1) vector containing weather each subject is treatment or control.
#' @param strata Numeric vector indexing the matched pairs and the reservoir. Entrees indexed with zero, are added to the reservoir. Entrees indexed with the same number are part of the same strata.
#' @param use_concordant_pairs Flag indicating whether to include the concordant pairs and reservoir in the model.
#' @param try_to_stableize Flag indicating whether the model should decompose the discordant pairs to increase stability if the model is unstable.
#' @return A new `clogitRClass` object.
#' @export
clogitR = function(response = NULL,
                   data = NULL,
                   treatment = NULL,
                   strata = NULL,
                   use_concordant_pairs = TRUE,
                   try_to_stableize = TRUE) {
  model = clogitRClass$new(
    response = response,
    data = data,
    treatment = treatment,
    strata = strata,
    use_concordant_pairs = use_concordant_pairs,
    try_to_stableize = try_to_stableize
  )
    
  class(model) = c("clogitR", class(model))
  model
}


#' @export
summary.clogitR = function(object, ...) {
  object$summary_data(...)  # just returns the summary object (no printing)
}

#' @export
print.summary.clogitR = function(x, ...) {
  
  
  cat(x$call, "\n\n")
  
  cat("Discordant model converged =", x$converged_discordant)
  if (x$use_concordant_pairs) { cat(" | Concordant and reservoir model converged =", x$converged_reservoir_concordant, '\n\n') } else {cat("\n\n")}
  
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
  
  # etc.
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
  
  
  
  
  