#' Initialize a new bclogit model
#'
#' @param response A binary (0,1) vector containing the response of each subject.
#' @param data A data.frame, data.table, or model.matrix containing the variables.
#' @param treatment A binary (0,1) vector containing whether each subject is treatment or control.
#' @param strata Numeric vector indexing the matched pairs and the reservoir. Entries indexed with zero, are added to the reservoir. Entries indexed with the same number are part of the same strata.
#' @param concordant_method The method to use for fitting the concordant pairs and reservoir. Options are "GLM", "GEE", and "GLMM".
#' @param prior_type The type of prior to use for the discordant pairs. Options are "naive", "G prior", "PMP", and "hybrid".
#' @return A list of class `bclogit`.
#' @export
bclogit <- function(x, ...) {
  UseMethod("bclogit")
}

#' @param treatment_name Optional string name for the treatment variable.
#' @describeIn bclogit Default method
#' @export
bclogit.default <- function(response = NULL,
                    data = NULL,
                    treatment = NULL,
                    strata = NULL,
                    concordant_method = "GEE",
                    prior_type = "hybrid",
                    treatment_name = NULL) {
  # ------------------------------------------------------------------------
  # 1. Input Validation and Pre-processing
  # ------------------------------------------------------------------------
  if (test_data_frame(data, types = c("numeric", "integer", "factor", "logical"))) {
    data_mat <- model.matrix(~ 0 + ., data = data)
  } else if (test_matrix(data, mode = "numeric")) {
    if (sd(data[, 1]) == 0) {
      data[, 1] <- NULL
    }
    data_mat <- as.matrix(data)
  } else {
    assert(
      check_data_frame(data),
      check_matrix(data, mode = "numeric")
    )
  }

  n <- nrow(data_mat)

  assertString(concordant_method, c("GLM", "GEE", "GLMM"))
  assertNumeric(response, lower = 0, upper = 1, any.missing = FALSE, len = n)

  if (!is.null(treatment)) {
    assertNumeric(treatment, lower = 0, upper = 1, any.missing = FALSE, len = n)
    if (!all(treatment %in% c(0, 1))) {
      stop("Treatment must be binary 0 or 1.")
    }
    if (is.null(treatment_name)) {
        treatment_name <- deparse(substitute(treatment))
    }
  }
  assertNumeric(strata, any.missing = FALSE, len = n)

  # Capture terms for summary/printing usage if possible (though we used matrix for fitting)
  model_terms <- tryCatch(terms(response ~ ., data = as.data.frame(data_mat)), error = function(e) NULL)

  # ------------------------------------------------------------------------
  # 2. Data Preparation
  # ------------------------------------------------------------------------

  X_concordant <- NULL
  y_concordant <- NULL
  treatment_concordant <- NULL

  X_diffs_discordant <- NULL
  y_diffs_discordant <- NULL
  treatment_diffs_discordant <- NULL

  # Internal variables for processing
  y <- response
  X <- data_mat
  X_model_matrix_col_names <- colnames(X)
  if (is.null(X_model_matrix_col_names)) {
    X_model_matrix_col_names <- paste0("X", 1:ncol(X))
  }

  if (length(y) <= ncol(X) + 5) {
    stop("Not enough rows. Must be at least the number of covariates plus 5")
  }

  if (is.null(treatment)) {
    warning("Function did not include a treatment parameter. The order of the differences in the matched pairs will be arbitrary.")
  }

  matched_data <- process_matched_pairs_cpp(
    strata = strata,
    y = y,
    X = X,
    treatment = treatment
  )

  X_concordant <- matched_data$X_concordant
  y_concordant <- matched_data$y_concordant
  treatment_concordant <- matched_data$treatment_concordant
  strata_concordant <- matched_data$strata_concordant
  X_diffs_discordant <- matched_data$X_diffs_discordant
  y_diffs_discordant <- matched_data$y_diffs_discordant
  treatment_diffs_discordant <- matched_data$treatment_diffs_discordant

  # ------------------------------------------------------------------------
  # 3. Model Fitting
  # ------------------------------------------------------------------------


  # --- Concordant Pairs / Reservoir Model ---

  concordant_model <- NULL
  concordant_converged <- NULL

  # Check for concordant pairs
  if (length(y_concordant) < ncol(X_concordant) + 5) {
    # Not enough data for concordant model
    warning("There are not enough concordant pairs or reservoir entries. The prior for the discordant pairs will be non-informative.")
    b_con <- rep(0, ncol(X_concordant) + 1)
    Sigma_con <- diag(101, ncol(X_concordant) + 1)
  } else {
    if (concordant_method == "GLM") {
      concordant_model <- glm(y_concordant ~ treatment_concordant + X_concordant, family = "binomial")
      b_con <- summary(concordant_model)$coefficients[, 1]
      Sigma_con <- pmin(vcov(concordant_model), 100)
      eps <- 1e-6 # Stability
      Sigma_con <- Sigma_con + diag(eps, nrow(Sigma_con))

      b_con <- c(0, b_con[-c(1, 2)])
      Sigma_con <- Sigma_con[-1, -1]
      Sigma_con[1, ] <- 0
      Sigma_con[, 1] <- 0
      Sigma_con[1, 1] <- 100
    } else if (concordant_method == "GEE") {
      concordant_model <- geepack::geeglm(
        y_concordant ~ treatment_concordant + X_concordant,
        id = strata_concordant,
        family = binomial(link = "logit"),
        corstr = "exchangeable",
        data = data.frame(y_concordant, treatment_concordant, X_concordant, strata_concordant)
      )

      b_con <- summary(concordant_model)$coefficients[, 1]
      Sigma_con <- pmin(vcov(concordant_model), 100)
      Sigma_con <- (Sigma_con + t(Sigma_con)) / 2 # Force symmetry
      eps <- 1e-6 # Stability
      Sigma_con <- Sigma_con + diag(eps, nrow(Sigma_con))

      b_con <- c(0, b_con[-c(1, 2)])
      Sigma_con <- Sigma_con[-1, -1]
      Sigma_con[1, ] <- 0
      Sigma_con[, 1] <- 0
      Sigma_con[1, 1] <- 100
    } else if (concordant_method == "GLMM") {
      concordant_model <- glmmTMB::glmmTMB(
        y_concordant ~ treatment_concordant + X_concordant + (1 | strata_concordant),
        family = binomial(),
        data   = data.frame(y_concordant, treatment_concordant, X_concordant, strata_concordant)
      )

      b_con <- summary(concordant_model)$coefficients$cond[, 1]
      Sigma_con <- pmin(vcov(concordant_model)$cond, 100)
      eps <- 1e-6 # Stability
      Sigma_con <- Sigma_con + diag(eps, nrow(Sigma_con))

      b_con <- c(0, b_con[-c(1, 2)])
      Sigma_con <- Sigma_con[-1, -1]
      Sigma_con[1, ] <- 0
      Sigma_con[, 1] <- 0
      Sigma_con[1, 1] <- 100
    }
  }


  discordant_model <- NULL
  converged_discordant <- NULL

  # --- Discordant Pairs Model ---
  if (length(y_diffs_discordant) < ncol(X_diffs_discordant) + 5) {
    stop("There are not enough discordant pairs. The model will not be fit .")
  } else {
    y_dis_0_1 <- ifelse(y_diffs_discordant == -1, 0, 1)
    wX_dis <- cbind(treatment_diffs_discordant, X_diffs_discordant)

    if (prior_type == "naive") {
      data_list <- list(
        N = nrow(wX_dis),
        K = ncol(wX_dis),
        X = wX_dis,
        y = y_dis_0_1,
        mu = b_con, # example prior mean
        Sigma = Sigma_con # example covariance (wide prior)
      )
      discordant_model <- tryCatch(
        {
          rstan::stan(file = system.file("stan/mvn_logistic.stan", package = "bclogit"), data = data_list, refresh = 0, chains = 1)
        },
        error = function(e) {
          warning(sprintf("stan_glm failed: %s", e$message))
          NULL
        }
      )
    } else if (prior_type == "G prior") {
      data_list <- list(
        N = nrow(wX_dis),
        K = ncol(wX_dis),
        X = wX_dis,
        y = y_dis_0_1,
        mu = b_con,
        Sigma = Sigma_con
      )
      discordant_model <- tryCatch(
        {
          rstan::stan(file = system.file("stan/mvn_logistic_gprior.stan", package = "bclogit"), data = data_list, refresh = 0, chains = 1)
        },
        error = function(e) {
          warning(sprintf("stan_glm failed: %s", e$message))
          NULL
        }
      )
    } else if (prior_type == "PMP") {
      # Orthogonalize inputs
      proj_matrix <- X_diffs_discordant %*% solve(t(X_diffs_discordant) %*% X_diffs_discordant) %*% t(X_diffs_discordant)
      w_dis_ortho <- treatment_diffs_discordant - proj_matrix %*% treatment_diffs_discordant

      data_list <- list(
        N = nrow(X_diffs_discordant),
        P = ncol(X_diffs_discordant),
        y = y_dis_0_1,
        xw = as.vector(w_dis_ortho),
        X = X_diffs_discordant,
        mu_A = as.array(b_con[-1]),
        Sigma_A = as.matrix(Sigma_con[-1, -1, drop = FALSE])
      )

      discordant_model <- tryCatch(
        {
          rstan::stan(file = system.file("stan/mvn_logistic_PMP.stan", package = "bclogit"), data = data_list, refresh = 0, chains = 1)
        },
        error = function(e) {
          warning(sprintf("stan_glm failed: %s", e$message))
          NULL
        }
      )
    } else if (prior_type == "Hybrid") {
      proj_matrix <- X_diffs_discordant %*% solve(t(X_diffs_discordant) %*% X_diffs_discordant) %*% t(X_diffs_discordant)
      w_dis_ortho <- treatment_diffs_discordant - proj_matrix %*% treatment_diffs_discordant


      # Prepare the data list for Stan
      data_list <- list(
        N = nrow(X_diffs_discordant),
        P = ncol(X_diffs_discordant),
        y = y_dis_0_1,
        xw = as.vector(w_dis_ortho),
        X = X_diffs_discordant,
        mu_A = as.array(b_con[-1]),
        Sigma_A = as.matrix(Sigma_con[-1, -1, drop = FALSE])
      )

      discordant_model <- tryCatch(
        {
          rstan::stan(file = system.file("stan/mvn_logistic_Hybrid.stan", package = "bclogit"), data = data_list, refresh = 0, chains = 1)
        },
        error = function(e) {
          warning(sprintf("Hybrid PMP failed: %s", e$message))
          NULL
        }
      )
    }
  }

  # ------------------------------------------------------------------------
  # 4. Result Construction
  # ------------------------------------------------------------------------

  # Extract coefficients from the Stan model if available
  coefficients <- NULL
  var_cov <- NULL

  if (!is.null(discordant_model)) {
    # Extract posterior samples

    sims <- rstan::extract(discordant_model)

    if (prior_type %in% c("naive", "G prior")) {
      beta_post <- sims$beta
      col_names_final <- c(treatment_name, X_model_matrix_col_names)
    } else {
      # PMP / Hybrid
      beta_w_post <- as.vector(sims$beta_w)
      beta_nuis_post <- sims$beta_nuis

      beta_post <- cbind(beta_w_post, beta_nuis_post)
      col_names_final <- c(treatment_name, X_model_matrix_col_names)
    }

    # Calculate posterior means
    coefficients <- colMeans(beta_post)
    names(coefficients) <- col_names_final

    # Calculate posterior covariance
    var_cov <- cov(beta_post)
    # Assign names
    rownames(var_cov) <- col_names_final
    colnames(var_cov) <- col_names_final
  }

  # Construct result object
  res <- list(
    coefficients = coefficients,
    var = var_cov,
    model = discordant_model,
    prior_info = list(
      mu = if (exists("b_con")) b_con else NULL,
      Sigma = if (exists("Sigma_con")) Sigma_con else NULL
    ),

    # Standard S3 components
    call = match.call(),
    terms = model_terms,
    xlevels = NULL,

    # Metadata for summary
    n = n,
    num_discordant = length(y_diffs_discordant),
    num_concordant = length(y_concordant),
    X_model_matrix_col_names = X_model_matrix_col_names,
    treatment_name = treatment_name
  )

  class(res) <- c("bclogit", "list")
  return(res)
}
