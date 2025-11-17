#include <RcppEigen.h>
using namespace Rcpp;

inline double tau(double z) {
  return 1.0 / (1.0 + std::exp(-z));
}

double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(Eigen::MatrixXd M, int j) {
  Eigen::VectorXd b;
  b.resize(M.rows());
  b.setZero();
  b(j - 1) = 1;

  Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower|Eigen::Upper> cg;
  cg.compute(M);
  Eigen::VectorXd x = cg.solve(b);

  return x(j - 1);
}

// Newton–Raphson solver on ΔX, y ∈ {±1}
// [[Rcpp::export(rng = false)]]
List fast_conditional_logistic_regression_cpp(const Eigen::MatrixXd &X_diff,
                                                        const Eigen::VectorXi &y_diff,
                                                        int max_iter = 100,
                                                        double tol = 1e-8) {
  int p = X_diff.cols();
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
  
  bool converged = false;
  for (int iter = 0; iter < max_iter; iter++) {
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(p);
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(p, p);

    for (int i = 0; i < X_diff.rows(); i++) {
      Eigen::VectorXd dx = X_diff.row(i);
      int outcome = y_diff(i);  // should be ±1

      double xb = dx.dot(beta);
      double p1 = tau(outcome * xb);
      double w  = p1 * (1.0 - p1);

      grad += (1.0 - p1) * outcome * dx;
      H    -= w * (dx * dx.transpose());
    }

    Eigen::VectorXd step = H.ldlt().solve(grad);
    beta -= step;
    if (step.norm() < tol){
      converged = true;
      break;
    }
  }

  return List::create(
      _["b"] = beta,
      _["converged"] = converged
    );
}

// Extended version: β + variance for coefficient j
// [[Rcpp::export(rng = false)]]
List fast_conditional_logistic_regression_with_var_cpp(const Eigen::MatrixXd &X_diff,
                                                      const Eigen::VectorXi &y_diff,
                                                      int j = 1,
                                                      int max_iter = 100,
                                                      double tol = 1e-8) {
  int p = X_diff.cols();
  List result = fast_conditional_logistic_regression_cpp(X_diff, y_diff, max_iter, tol);
  Eigen::VectorXd beta = as<Eigen::VectorXd>(result["b"]);
  // Fisher info at β̂
  Eigen::MatrixXd I = Eigen::MatrixXd::Zero(p, p);
  for (int i = 0; i < X_diff.rows(); i++) {
    Eigen::VectorXd dx = X_diff.row(i);
    int outcome = y_diff(i);

    double xb = dx.dot(beta);
    double p1 = tau(outcome * xb);
    double w  = p1 * (1.0 - p1);

    // ADD SMALL REGULARIZATION TO PREVENT COMPLETE SEPARATION ISSUES
    w = std::max(w, 1e-10); 

    I += w * (dx * dx.transpose());
  }

  // ADD RIDGE REGULARIZATION TO FISHER INFORMATION MATRIX
  double ridge = 1e-6;
  I += ridge * Eigen::MatrixXd::Identity(p, p);

  // Use LDLT decomposition instead of direct inversion - much faster and more stable
  Eigen::LDLT<Eigen::MatrixXd> ldlt(I);
  
  // Check if decomposition succeeded
  if (!ldlt.isPositive()) {
    // Matrix is singular or not positive definite - return NaN
    if (j == -1){
      Eigen::VectorXd nan_vec = Eigen::VectorXd::Constant(p, std::numeric_limits<double>::quiet_NaN());
      return List::create(
        _["b"] = beta,
        _["ssq_b"] = nan_vec,
        _["converged"] = result["converged"]
      );
    } else if (j == -2){
      Eigen::MatrixXd nan_mat = Eigen::MatrixXd::Constant(p, p, std::numeric_limits<double>::quiet_NaN());
      return List::create(
        _["b"] = beta,
        _["ssq_b"] = nan_mat,
        _["converged"] = result["converged"]
      );
    } else {
      return List::create(
        _["b"] = beta,
        _["ssq_b_j"] = std::numeric_limits<double>::quiet_NaN(),
        _["converged"] = result["converged"]
      );
    }
  }

  if (j == -1){
    // Compute diagonal of inverse using LDLT solve
    Eigen::VectorXd diag_inv(p);
    for (int i = 0; i < p; i++) {
      Eigen::VectorXd e = Eigen::VectorXd::Zero(p);
      e(i) = 1.0;
      diag_inv(i) = ldlt.solve(e)(i);
    }
    return List::create(
      _["b"] = beta,
      _["ssq_b"] = diag_inv,
      _["converged"] = result["converged"]
    );
  } else if (j == -2){
    // Compute full inverse using LDLT solve
    Eigen::MatrixXd inv = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
    return List::create(
      _["b"] = beta,
      _["ssq_b"] = inv,
      _["converged"] = result["converged"]
    );
  } else {
    // Compute single diagonal entry
    Eigen::VectorXd e = Eigen::VectorXd::Zero(p);
    e(j - 1) = 1.0;
    double ssq_j = ldlt.solve(e)(j - 1);
    return List::create(
      _["b"] = beta,
      _["ssq_b_j"] = ssq_j,
      _["converged"] = result["converged"]
    );
  }
}
