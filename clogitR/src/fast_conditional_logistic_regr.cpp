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
// [[Rcpp::export]]
Eigen::VectorXd fast_conditional_logistic_regression_cpp(const Eigen::MatrixXd &X_diff,
                                                        const Eigen::VectorXi &y_diff,
                                                        int max_iter = 100,
                                                        double tol = 1e-8) {
  int p = X_diff.cols();
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);

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
    if (step.norm() < tol) break;
  }

  return beta;
}

// Extended version: β + variance for coefficient j
// [[Rcpp::export]]
List fast_conditional_logistic_regression_with_var_cpp(const Eigen::MatrixXd &X_diff,
                                                      const Eigen::VectorXi &y_diff,
                                                      int j = 1,
                                                      int max_iter = 100,
                                                      double tol = 1e-8) {
  int p = X_diff.cols();
  Eigen::VectorXd beta = fast_conditional_logistic_regression_cpp(X_diff, y_diff, max_iter, tol);

  // Fisher info at β̂
  Eigen::MatrixXd I = Eigen::MatrixXd::Zero(p, p);
  for (int i = 0; i < X_diff.rows(); i++) {
    Eigen::VectorXd dx = X_diff.row(i);
    int outcome = y_diff(i);

    double xb = dx.dot(beta);
    double p1 = tau(outcome * xb);
    double w  = p1 * (1.0 - p1);

    I += w * (dx * dx.transpose());
  }

  if (j == -1){
    return List::create(
      _["b"] = beta,
      _["ssq_b"] = I.solve().diagonal()
    );
  } else {
    return List::create(
      _["b"] = beta,
      _["ssq_b_j"] = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(I, j)
    );
  }
}
