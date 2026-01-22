
# Validating S3 Methods: coef, vcov, confint, formula
devtools::load_all("bclogit")

# Generate dummy matched pairs data
set.seed(456)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2)
colnames(X) <- c("age", "weight")
treatment <- rep(c(0, 1), n / 2)
strata <- rep(1:(n / 2), each = 2)
# Linear predictor
eta <- -0.5 + 1.5 * treatment + 0.5 * X[, 1] - 0.5 * X[, 2]
probs <- 1 / (1 + exp(-eta))
response <- rbinom(n, 1, probs)

cat("\n--- Fitting Model ---\n")
# Use naive prior for speed in test, and GLM for stability on random data
fit <- bclogit(response = response, data = X, treatment = treatment, strata = strata, 
              prior_type = "naive", concordant_method = "GLM")

cat("\n--- Testing coef() ---\n")
c <- coef(fit)
print(c)
if (is.null(c)) stop("coef() returned NULL")
if (length(c) != 3) stop("coef() wrong length (expected 3: treatment, age, weight)")

cat("\n--- Testing vcov() ---\n")
v <- vcov(fit)
print(v)
if (!is.matrix(v) || nrow(v) != 3 || ncol(v) != 3) stop("vcov() returned invalid matrix")

cat("\n--- Testing formula() ---\n")
f <- formula(fit)
print(f)
if (class(f) != "formula") stop("formula() did not return a formula object")

cat("\n--- Testing confint() ---\n")
ci <- confint(fit, level = 0.95)
print(ci)
if (nrow(ci) != 3 || ncol(ci) != 2) stop("confint() returned invalid dimensions")
if (any(is.na(ci))) stop("confint() returned NA values")

cat("\n--- Testing confint() subset ---\n")
ci_sub <- confint(fit, parm = "treatment", level = 0.90)
print(ci_sub)
if (nrow(ci_sub) != 1) stop("confint() subset returned wrong rows")

cat("\nS3 Methods Test Passed\n")
