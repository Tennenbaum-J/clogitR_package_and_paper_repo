
#' @export
#' @describeIn bclogit Formula method
bclogit.formula <- function(formula, data, treatment = NULL, strata = NULL, ...) {
  cl <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "drop.unused.levels"), names(cl), 0L)
  cl <- cl[c(1L, m)]
  cl[[1L]] <- quote(stats::model.frame)
  mf <- eval(cl, parent.frame())
  
  response <- model.response(mf, "numeric")
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)
  
  if (colnames(X)[1] == "(Intercept)") {
      X <- X[, -1, drop = FALSE]
  }
  
  if (missing(data)) data <- environment(formula)

  trt_vec <- eval(substitute(treatment), data, parent.frame())
  strata_vec <- eval(substitute(strata), data, parent.frame())
  
  # Extract treatment name from call
  # We use the full match.call() which has all arguments
  full_cl <- match.call()
  if ("treatment" %in% names(full_cl)) {
      t_name <- deparse(full_cl$treatment)
      # Clean up if it's long or complex? usually deparse is fine.
  } else {
      t_name <- "treatment"
  }
  
  bclogit.default(response = response, data = X, treatment = trt_vec, strata = strata_vec, treatment_name = t_name, ...)
}

#' @export
summary.bclogit <- function(object, conf.level = 0.95, ...) {
  if (is.null(object$model)) {
    # Fallback if no model was fit / converged
    cat("No discordant model available.\n")
    return(invisible(object))
  }

  # Extract posterior samples again or use stored summary if optimization needed
  # We'll extract from standard 'stanfit' object
  sims <- rstan::extract(object$model)

  # Reconstruct beta_post matrix
  if (!is.null(sims$beta)) {
    beta_post <- sims$beta
  } else {
    beta_post <- cbind(as.vector(sims$beta_w), sims$beta_nuis)
  }

  # Calculate summary stats
  est <- colMeans(beta_post)
  se <- apply(beta_post, 2, sd)

  alpha <- (1 - conf.level) / 2
  q_low <- apply(beta_post, 2, quantile, probs = alpha)
  q_high <- apply(beta_post, 2, quantile, probs = 1 - alpha)

  # Probability of being positive/negative (for hypothesis testing)
  # Pr(beta > 0)
  prob_pos <- colMeans(beta_post > 0)

  coef_mat <- cbind(
    Estimate = est,
    `Est.Error` = se,
    `Q_low` = q_low,
    `Q_high` = q_high,
    `Pr(>0)` = prob_pos
  )

  # Set row names
  rownames(coef_mat) <- names(object$coefficients)

  # Rename Quantile columns nicely
  colnames(coef_mat)[3:4] <- paste0(c("L", "U"), format(conf.level * 100), "%")

  res <- list(
    call = object$call,
    coefficients = coef_mat,
    num_discordant = object$num_discordant,
    num_concordant = object$num_concordant,
    conf.level = conf.level
  )

  class(res) <- "summary.bclogit"
  res
}

#' @export
print.summary.bclogit <- function(x, digits = 4, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\n")

  cat("Discordant Pairs Bayesian Model\n")
  cat("Number of discordant pairs:", x$num_discordant, "\n")
  cat("Number of concordant or reservoir entries (used for prior):", x$num_concordant, "\n")
  cat("\n")

  cat("Coefficients (Posterior Mean and Quantiles):\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = FALSE)

  cat("\nNote: Pr(>0) is the posterior probability that the coefficient is positive.\n")

  invisible(x)
}

#' @export
coef.bclogit <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.bclogit <- function(object, ...) {
  object$var
}

#' @export
formula.bclogit <- function(x, ...) {
  if (!is.null(x$terms)) {
    formula(x$terms)
  } else {
    NULL
  }
}

#' @export
confint.bclogit <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$model)) {
    warning("No discordant model available for confidence intervals.")
    return(NULL)
  }

  sims <- rstan::extract(object$model)
  if (!is.null(sims$beta)) {
    beta_post <- sims$beta
  } else {
    beta_post <- cbind(as.vector(sims$beta_w), sims$beta_nuis)
  }
  
  if (ncol(beta_post) == length(object$coefficients)) {
      colnames(beta_post) <- names(object$coefficients)
  }

  if (missing(parm)) {
    parm <- names(object$coefficients)
  } else if (is.numeric(parm)) {
    parm <- names(object$coefficients)[parm]
  }
  
  present_parms <- intersect(parm, colnames(beta_post))
  if (length(present_parms) == 0) {
      stop("No parameters found matching the 'parm' argument.")
  }

  alpha <- (1 - level) / 2
  probs <- c(alpha, 1 - alpha)
  
  ci <- apply(beta_post[, present_parms, drop = FALSE], 2, quantile, probs = probs)
  t(ci)
}
