#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

class Evaluator {
private:
  List true_data;
  List imputed_data;
  std::string method;

  List metrics; // column -> vector of metrics
  double rmse, mae, r2, corr, ks, acc; // aggregated global metrics

  // Helper: compute metrics for one column
  NumericVector compute_column_metrics(NumericVector true_vals, NumericVector imputed_vals) {
    int n = true_vals.size();
    if (n != imputed_vals.size()) {
      stop("true and imputed vectors must have same length");
    }

    // RMSE & MAE
    double se = 0.0, abs_err = 0.0;
    double ss_res = 0.0, ss_tot = 0.0;
    double mean_true = mean(true_vals);
    int correct_count = 0;

    for (int i = 0; i < n; i++) {
      double diff = true_vals[i] - imputed_vals[i];
      se += diff * diff;
      abs_err += std::abs(diff);
      ss_res += diff * diff;
      ss_tot += (true_vals[i] - mean_true) * (true_vals[i] - mean_true);

      if (true_vals[i] == imputed_vals[i]) correct_count++;
    }

    double col_rmse = std::sqrt(se / n);
    double col_mae = abs_err / n;
    double col_r2 = (ss_tot == 0) ? NA_REAL : 1.0 - (ss_res / ss_tot);

    // Fixed correlation calculation
    double col_corr = NA_REAL;
    {
      double mean_x = mean(true_vals);
      double mean_y = mean(imputed_vals);
      double num = 0.0, den_x = 0.0, den_y = 0.0;
      for (int i = 0; i < n; i++) {
        double dx = true_vals[i] - mean_x;
        double dy = imputed_vals[i] - mean_y;
        num += dx * dy;
        den_x += dx * dx;
        den_y += dy * dy;
      }
      if (den_x > 0 && den_y > 0) {
        col_corr = num / std::sqrt(den_x * den_y);
      }
    }

    // KS statistic
    double col_ks = ks_test(true_vals, imputed_vals);

    // Accuracy (exact match proportion)
    double col_acc = (double) correct_count / n;

    NumericVector results = NumericVector::create(
      _["RMSE"] = col_rmse,
      _["MAE"] = col_mae,
      _["R2"] = col_r2,
      _["Correlation"] = col_corr,
      _["KS"] = col_ks,
      _["Accuracy"] = col_acc
    );

    return results;
  }

  // Kolmogorov-Smirnov statistic (simple implementation)
  double ks_test(NumericVector x, NumericVector y) {
    NumericVector xs = clone(x).sort();
    NumericVector ys = clone(y).sort();

    int nx = xs.size(), ny = ys.size();
    int i = 0, j = 0;
    double cdf_x = 0.0, cdf_y = 0.0, d = 0.0;

    while (i < nx && j < ny) {
      if (xs[i] <= ys[j]) {
        i++;
        cdf_x = (double) i / nx;
      } else {
        j++;
        cdf_y = (double) j / ny;
      }
      d = std::max(d, std::abs(cdf_x - cdf_y));
    }
    return d;
  }

  void compute_all_metrics() {
    CharacterVector cols = true_data.names();
    int p = cols.size();

    NumericVector rmse_vals(p), mae_vals(p), r2_vals(p), corr_vals(p), ks_vals(p), acc_vals(p);

    for (int i = 0; i < p; i++) {
      std::string cname = as<std::string>(cols[i]);
      NumericVector tcol = true_data[cname];
      NumericVector icol = imputed_data[cname];

      NumericVector res = compute_column_metrics(tcol, icol);

      metrics[cname] = res;

      rmse_vals[i] = res["RMSE"];
      mae_vals[i] = res["MAE"];
      r2_vals[i] = res["R2"];
      corr_vals[i] = res["Correlation"];
      ks_vals[i] = res["KS"];
      acc_vals[i] = res["Accuracy"];
    }

    rmse = mean(rmse_vals);
    mae = mean(mae_vals);
    r2 = mean(r2_vals);
    corr = mean(corr_vals);
    ks = mean(ks_vals);
    acc = mean(acc_vals);
  }

public:
  Evaluator(List true_data_, List imputed_data_, std::string method_)
    : true_data(true_data_), imputed_data(imputed_data_), method(method_) {
    compute_all_metrics();
  }

  List toList() {
    return List::create(
      _["method"] = method,
      _["metrics"] = metrics, // per-column dictionary
      _["RMSE"] = rmse,
      _["MAE"] = mae,
      _["R2"] = r2,
      _["Correlation"] = corr,
      _["KS"] = ks,
      _["Accuracy"] = acc
    );
  }
};

// [[Rcpp::export]]
List evaluate_imputation(List true_data, List imputed_data, std::string method) {
  Evaluator eval(true_data, imputed_data, method);
  return eval.toList();
}
