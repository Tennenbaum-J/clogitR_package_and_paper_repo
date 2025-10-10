#' Initialize a new clogitR model
#'
#' @export
#' @param formula A formula specifying the model.
#' @param data A data.frame or data.table containing the variables.
#' @param use_concordant_pairs Flag indicating whether to include the concordant pairs and reservoir in the model.
#' @param try_to_stableize Flag indicating whether the model should decompose the discordant pairs to increase stability if the model is unstable.
#' @return A new `clogitRClass` object.
clogitR = function(formula, data,
                    use_concordant_pairs = TRUE,
                    try_to_stableize = TRUE) {
  model = clogitRClass$new(
    formula = formula,
    data = data,
    use_concordant_pairs = use_concordant_pairs,
    try_to_stableize = try_to_stableize
  )
  class(model) = c("clogitR", class(model))
  model
}


#' @export
summary.clogitR = function(object, ...) {
  object$summary(...)
}

