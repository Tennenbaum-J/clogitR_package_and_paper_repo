#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export(rng = false)]]
List process_matched_pairs_cpp(
    const Eigen::VectorXi& strata,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& X,
    Nullable<Eigen::VectorXd> treatment = R_NilValue
) {
  const int n = strata.size();
  const int p = X.cols();
  const bool has_treatment = treatment.isNotNull();
  
  Eigen::VectorXd treat_vec;
  if (has_treatment) {
    treat_vec = as<Eigen::VectorXd>(treatment);
  }
  
  // Pre-count to pre-allocate (single pass through strata)
  int n_reservoir = 0;
  int n_discordant = 0;
  int max_strata = 0;
  
  for (int i = 0; i < n; i++) {
    if (strata[i] == 0) {
      n_reservoir++;
    }
    if (strata[i] > max_strata) {
      max_strata = strata[i];
    }
  }
  
  // Pre-allocate with upper bounds (concordant pairs will add to reservoir)
  int max_concordant = max_strata; // worst case: all pairs are concordant
  Eigen::MatrixXd reservoir_X(n_reservoir + 2 * max_concordant, p);
  Eigen::VectorXd reservoir_y(n_reservoir + 2 * max_concordant);
  Eigen::VectorXd reservoir_treat(has_treatment ? n_reservoir + 2 * max_concordant : 0);
  
  Eigen::MatrixXd diffs_X(max_strata, p); // at most max_strata discordant pairs
  Eigen::VectorXd diffs_y(max_strata);
  Eigen::VectorXd diffs_treat(has_treatment ? max_strata : 0);
  std::vector<int> discordant_idx;
  
  
  int res_idx = 0;
  int diff_idx = 0;
  
  // Handle reservoir (strata == 0) - fast vectorized operation
  for (int i = 0; i < n; i++) {
    if (strata[i] == 0) {
      reservoir_X.row(res_idx) = X.row(i);
      reservoir_y[res_idx] = y[i];
      if (has_treatment) {
        reservoir_treat[res_idx] = treat_vec[i];
      }
      res_idx++;
    }
  }
  
  // Build index for fast pair lookup - O(n) instead of O(n*max_strata)
  std::vector<std::vector<int>> pair_indices(max_strata + 1);
  for (int i = 0; i < n; i++) {
    if (strata[i] > 0) {
      pair_indices[strata[i]].push_back(i);
    }
  }
  
  // Process all pairs
  for (int pair_num = 1; pair_num <= max_strata; pair_num++) {
    const std::vector<int>& pair = pair_indices[pair_num];
    
    if (pair.size() == 0) {
      stop("Stratum index numbers must be sequential.");
    }
    if (pair.size() != 2) {
      stop("Each nonzero stratum must have exactly 2 rows.");
    }
    
    const int i = pair[0];
    const int j = pair[1];
    const double yi = y[i];
    const double yj = y[j];
    
    if (yi == yj) {
      // Concordant pair → add to reservoir
      reservoir_X.row(res_idx) = X.row(i);
      reservoir_y[res_idx] = yi;
      if (has_treatment) {
        reservoir_treat[res_idx] = treat_vec[i];
      }
      res_idx++;
      
      reservoir_X.row(res_idx) = X.row(j);
      reservoir_y[res_idx] = yj;
      if (has_treatment) {
        reservoir_treat[res_idx] = treat_vec[j];
      }
      res_idx++;
    } else {
      // Discordant → compute diff (case - control)
      if (has_treatment && treat_vec[i] == 0) {
        discordant_idx.push_back(j);
        discordant_idx.push_back(i);
        diffs_X.row(diff_idx) = X.row(j) - X.row(i);
        diffs_y[diff_idx] = yj - yi;
        diffs_treat[diff_idx] = treat_vec[j] - treat_vec[i];
      } else {
        discordant_idx.push_back(i);
        discordant_idx.push_back(j);
        diffs_X.row(diff_idx) = X.row(i) - X.row(j);
        diffs_y[diff_idx] = yi - yj;
        if (has_treatment) {
          diffs_treat[diff_idx] = treat_vec[i] - treat_vec[j];
        }
      }
      diff_idx++;
    }
  }
  
  // Check if we have enough discordant pairs
  if (diff_idx < p + 5) {
    warning("There are not enough discordant pairs. All the data was put in the reservoir.");
    return List::create(
      _["X_reservoir_concordant"] = X,
      _["y_reservoir_concordant"] = y,
      _["treatment_reservoir_concordant"] = treatment,
      _["X_diffs_discordant"] = R_NilValue,
      _["y_diffs_discordant"] = R_NilValue,
      _["treatment_diffs_discordant"] = R_NilValue,
      _["dropped_discordant"] = true,
      _["dropped_reservoir_concordant"] = false,
      _["discordant_idx"] = R_NilValue
    );
  }
  
  bool dropped_reservoir = (res_idx <= p + 5);
  
  // Resize to actual sizes (cheap operation, just changes dimensions)
  if (res_idx > 0) {
    reservoir_X.conservativeResize(res_idx, p);
    reservoir_y.conservativeResize(res_idx);
    if (has_treatment) {
      reservoir_treat.conservativeResize(res_idx);
    }
  }
  
  if (diff_idx > 0) {
    diffs_X.conservativeResize(diff_idx, p);
    diffs_y.conservativeResize(diff_idx);
    if (has_treatment) {
      diffs_treat.conservativeResize(diff_idx);
    }
  }

  return List::create(
    _["X_reservoir_concordant"] = res_idx > 0 ? wrap(reservoir_X) : R_NilValue,
    _["y_reservoir_concordant"] = res_idx > 0 ? wrap(reservoir_y) : R_NilValue,
    _["treatment_reservoir_concordant"] = (has_treatment && res_idx > 0) ? wrap(reservoir_treat) : R_NilValue,
    _["X_diffs_discordant"] = diff_idx > 0 ? wrap(diffs_X) : R_NilValue,
    _["y_diffs_discordant"] = diff_idx > 0 ? wrap(diffs_y) : R_NilValue,
    _["treatment_diffs_discordant"] = (has_treatment && diff_idx > 0) ? wrap(diffs_treat) : R_NilValue,
    _["dropped_discordant"] = false,
    _["dropped_reservoir_concordant"] = dropped_reservoir,
    _["discordant_idx"] = discordant_idx
  );
}