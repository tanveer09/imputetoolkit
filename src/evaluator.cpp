#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

class Evaluator {
private:
  NumericVector true_vals;
  NumericVector imputed_vals;
  double rmse;
  double mae;
  double r2;
  std::string method;

  void compute_metrics() {
    int n = true_vals.size();
    if (n != imputed_vals.size()) {
      stop("true_vals and imputed_vals must have the same length");
    }

    double ss_res = 0.0, ss_tot = 0.0, mean_true = mean(true_vals);
    double abs_error_sum = 0.0;
    double se = 0.0;

    for (int i = 0; i < n; i++) {
      double diff = true_vals[i] - imputed_vals[i];
      se += diff * diff;
      abs_error_sum += std::abs(diff);
      ss_res += diff * diff;
      ss_tot += (true_vals[i] - mean_true) * (true_vals[i] - mean_true);
    }

    rmse = std::sqrt(se / n);
    mae = abs_error_sum / n;
    r2 = 1.0 - (ss_res / ss_tot);
  }

public:
  Evaluator(NumericVector true_vals_, NumericVector imputed_vals_, std::string method_)
    : true_vals(true_vals_), imputed_vals(imputed_vals_), method(method_) {
    compute_metrics();
  }

  List toList() {
    return List::create(
      _["method"] = method,
      _["RMSE"]   = rmse,
      _["MAE"]    = mae,
      _["R2"]     = r2
    );
  }
};

// [[Rcpp::export]]
List evaluate_imputation(NumericVector true_vals, NumericVector imputed_vals, std::string method) {
  Evaluator eval(true_vals, imputed_vals, method);
  return eval.toList();
}
