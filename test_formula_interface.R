
# Validation: Formula Interface
devtools::load_all("bclogit")

set.seed(123)
n <- 200
df <- data.frame(
  age = rnorm(n),
  weight = rnorm(n),
  treatment = rep(c(0, 1), n / 2),
  strata = rep(1:(n / 2), each = 2)
)

eta <- -0.5 + 1.2 * df$treatment + 0.4 * df$age - 0.3 * df$weight
probs <- 1 / (1 + exp(-eta))
df$y <- rbinom(n, 1, probs)

cat("\n--- Fitting Formula Model ---\n")
# bclogit(formula, data, treatment, strata, ...)
# Note: treatment and strata must be columns in data or vectors in environment
fit <- bclogit(y ~ age + weight, data = df, treatment = treatment, strata = strata, 
               prior_type = "naive", concordant_method = "GLM")

print(fit)
s <- summary(fit)
print(s)

cat("\n--- Checking Coefficients ---\n")
c <- coef(fit)
print(c)
if (!all(c("treatment", "age", "weight") %in% names(c))) {
    stop("Coefficients missing expected names.")
}

cat("\n--- Checking Formula Object ---\n")
f <- formula(fit)
print(f)
# Should be y ~ age + weight
if (format(f) != "y ~ age + weight") {
   warning(paste("Formula mismatch:", format(f)))
}

cat("\nFormula Interface Test Passed\n")
