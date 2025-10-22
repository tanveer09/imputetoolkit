// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <map>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// -------------------------
// --- Helper: KS test ---
// -------------------------
double ks_test(NumericVector x, NumericVector y) {
  NumericVector xs = clone(x).sort();
  NumericVector ys = clone(y).sort();
  int nx = xs.size(), ny = ys.size();
  int i = 0, j = 0;
  double cdf_x = 0.0, cdf_y = 0.0, d = 0.0;

  while (i < nx && j < ny) {
    if (xs[i] <= ys[j]) {
      i++; cdf_x = (double)i / nx;
    } else {
      j++; cdf_y = (double)j / ny;
    }
    d = std::max(d, std::abs(cdf_x - cdf_y));
  }
  return d;
}

// -------------------------
// --- Numeric metrics ---
// -------------------------
List compute_numeric_metrics(NumericVector y, NumericVector yhat) {
  if (y.size() != yhat.size())
    stop("Length mismatch in numeric metric computation");

  // --- Remove missing pairs ---
  LogicalVector ok = !is_na(y) & !is_na(yhat);
  NumericVector y_clean = y[ok];
  NumericVector yhat_clean = yhat[ok];
  int n = y_clean.size();

  // --- Handle empty or constant vectors safely ---
  if (n == 0) {
    return List::create(
      _["RMSE"] = 0.0,
      _["MAE"]  = 0.0,
      _["R2"]   = NA_REAL,
      _["Correlation"] = NA_REAL,
      _["KS"] = 0.0
    );
  }

  double se = 0.0, ae = 0.0;
  double ss_res = 0.0, ss_tot = 0.0;
  double mean_y = mean(y_clean);

  for (int i = 0; i < n; i++) {
    double diff = y_clean[i] - yhat_clean[i];
    se += diff * diff;
    ae += std::abs(diff);
    ss_res += diff * diff;
    ss_tot += (y_clean[i] - mean_y) * (y_clean[i] - mean_y);
  }

  double rmse = (n > 0) ? std::sqrt(se / n) : 0.0;
  double mae  = (n > 0) ? ae / n : 0.0;
  double r2   = (ss_tot == 0) ? NA_REAL : 1.0 - (ss_res / ss_tot);

  // --- Safe correlation ---
  // --- Safe correlation computation ---
  double corr = 0.0;  // default to 0 (no correlation)
  if (n > 1) {
    double mean_x = mean(y_clean);
    double mean_yhat = mean(yhat_clean);
    double num = 0.0, denx = 0.0, deny = 0.0;

    for (int i = 0; i < n; i++) {
      double dx = y_clean[i] - mean_x;
      double dy = yhat_clean[i] - mean_yhat;
      num  += dx * dy;
      denx += dx * dx;
      deny += dy * dy;
    }

    if (denx > 0 && deny > 0) {
      corr = num / std::sqrt(denx * deny);
      // Clamp correlation to [-1, 1] range in case of rounding noise
      if (!R_finite(corr)) corr = 0.0;
      if (corr > 1.0) corr = 1.0;
      if (corr < -1.0) corr = -1.0;
    } else {
      // one or both variables constant
      corr = 0.0;
    }
  } else {
    corr = 0.0;
  }


  double ks = (n > 1) ? ks_test(y_clean, yhat_clean) : 0.0;

  // --- Replace any remaining NAs with safe defaults ---
  if (NumericVector::is_na(rmse)) rmse = 0.0;
  if (NumericVector::is_na(mae)) mae = 0.0;
  if (NumericVector::is_na(ks)) ks = 0.0;

  return List::create(
    _["RMSE"] = rmse,
    _["MAE"]  = mae,
    _["R2"]   = r2,
    _["Correlation"] = corr,
    _["KS"] = ks
  );
}


// -------------------------
// --- Categorical metrics ---
// -------------------------
Rcpp::List compute_categorical_metrics(Rcpp::IntegerVector y, Rcpp::IntegerVector yhat) {
  int n = y.size();
  if (n != yhat.size()) Rcpp::stop("Length mismatch in categorical metric computation");

  int correct = 0;
  std::map<int,int> true_count, pred_count, match_count;

  for (int i = 0; i < n; i++) {
    true_count[y[i]]++;
    pred_count[yhat[i]]++;
    if (y[i] == yhat[i]) {
      correct++;
      match_count[y[i]]++;
    }
  }

  double accuracy = static_cast<double>(correct) / n;

  // Cohen's Kappa
  double pe = 0.0;
  for (auto &p : true_count) {
    int c = p.first;
    pe += (static_cast<double>(true_count[c]) / n) * (static_cast<double>(pred_count[c]) / n);
  }
  double kappa = ((1.0 - pe) == 0.0) ? NA_REAL : (accuracy - pe) / (1.0 - pe);

  // Per-class metrics for Macro-F1 and Balanced Accuracy
  double macro_f1 = 0.0;
  double macro_recall = 0.0;
  double classes = static_cast<double>(true_count.size());

  for (auto &p : true_count) {
    int c = p.first;
    double tp = match_count[c];
    double fp = pred_count[c] - tp;
    double fn = true_count[c] - tp;

    double precision = ((tp + fp) == 0.0) ? 0.0 : tp / (tp + fp);
    double recall    = ((tp + fn) == 0.0) ? 0.0 : tp / (tp + fn);
    double f1        = ((precision + recall) == 0.0) ? 0.0 : (2.0 * precision * recall) / (precision + recall);

    macro_f1     += f1;
    macro_recall += recall;
  }
  macro_f1     = (classes > 0.0) ? (macro_f1 / classes)     : NA_REAL;
  double balanced_acc = (classes > 0.0) ? (macro_recall / classes) : NA_REAL;

  // Micro-F1 equals Accuracy in single-label multi-class classification
  double micro_f1 = accuracy;


  return Rcpp::List::create(
    Rcpp::_["Accuracy"]      = accuracy,
    Rcpp::_["Kappa"]         = kappa,
    Rcpp::_["F1"]            = micro_f1,       // NEW (micro-F1)
    Rcpp::_["MacroF1"]       = macro_f1,
    Rcpp::_["BalancedAccuracy"]   = balanced_acc   // renamed to match test
  );
}


// -------------------------
// --- Main entry point ---
// -------------------------

// [[Rcpp::export]]
Rcpp::List evaluate_imputation_split(Rcpp::List true_data, Rcpp::List imputed_data,
                                     Rcpp::CharacterVector numeric_cols,
                                     Rcpp::CharacterVector categorical_cols,
                                     std::string method) {

  Rcpp::List numeric_metrics, categorical_metrics;

  // --- Numeric ---
  for (int i = 0; i < numeric_cols.size(); i++) {
    std::string colname = Rcpp::as<std::string>(numeric_cols[i]);

    // Skip if this column was not evaluated (no injected indices)
    if (!true_data.containsElementNamed(colname.c_str())) continue;
    if (!imputed_data.containsElementNamed(colname.c_str())) continue;

    Rcpp::NumericVector y    = true_data[colname];
    Rcpp::NumericVector yhat = imputed_data[colname];
    numeric_metrics[colname] = compute_numeric_metrics(y, yhat);
  }

  // --- Categorical ---
  for (int i = 0; i < categorical_cols.size(); i++) {
    std::string colname = Rcpp::as<std::string>(categorical_cols[i]);

    // Skip if not evaluated
    if (!true_data.containsElementNamed(colname.c_str())) continue;
    if (!imputed_data.containsElementNamed(colname.c_str())) continue;

    // Factors come in as integers in Rcpp; coerce safely
    Rcpp::IntegerVector y    = Rcpp::as<Rcpp::IntegerVector>(true_data[colname]);
    Rcpp::IntegerVector yhat = Rcpp::as<Rcpp::IntegerVector>(imputed_data[colname]);
    categorical_metrics[colname] = compute_categorical_metrics(y, yhat);
  }

  return Rcpp::List::create(
    Rcpp::_["method"]               = method,
    Rcpp::_["metrics_numeric"]      = numeric_metrics,
    Rcpp::_["metrics_categorical"]  = categorical_metrics
  );
}


// [[Rcpp::export]]
Rcpp::List evaluate_imputation(Rcpp::List true_data, Rcpp::List imputed_data, std::string method) {
  // Fallback simple evaluator for numeric-only data
  Rcpp::List res_numeric;

  for (int i = 0; i < true_data.size(); i++) {
    Rcpp::NumericVector t = true_data[i];
    Rcpp::NumericVector p = imputed_data[i];

    int n = t.size();
    if (n != p.size()) stop("Length mismatch in evaluate_imputation");

    double se = 0.0, ae = 0.0, ss_res = 0.0, ss_tot = 0.0;
    double mean_t = Rcpp::mean(t);

    for (int j = 0; j < n; j++) {
      double diff = t[j] - p[j];
      se += diff * diff;
      ae += std::abs(diff);
      ss_res += diff * diff;
      ss_tot += (t[j] - mean_t) * (t[j] - mean_t);
    }

    double rmse = std::sqrt(se / n);
    double mae  = ae / n;
    double r2   = (ss_tot == 0) ? NA_REAL : 1.0 - (ss_res / ss_tot);

    // Compute correlation manually (safe)
    double mean_p = Rcpp::mean(p);
    double num = 0.0, denx = 0.0, deny = 0.0;
    for (int j = 0; j < n; j++) {
      double dx = t[j] - mean_t;
      double dy = p[j] - mean_p;
      num  += dx * dy;
      denx += dx * dx;
      deny += dy * dy;
    }
    double correlation = (denx > 0 && deny > 0) ? num / std::sqrt(denx * deny) : NA_REAL;

    // KS statistic
    double ks = ks_test(t, p);

    // Accuracy for numeric equality (use tolerance)
    int equal = 0;
    for (int j = 0; j < n; j++) {
      if (std::abs(t[j] - p[j]) < 1e-8) equal++;
    }
    double accuracy = (double)equal / n;

    Rcpp::List m = Rcpp::List::create(
      _["RMSE"] = rmse,
      _["MAE"]  = mae,
      _["R2"]   = r2,
      _["Correlation"] = correlation,
      _["KS"] = ks,
      _["Accuracy"] = accuracy
    );

    res_numeric.push_back(m);
  }

  return Rcpp::List::create(
    _["method"] = method,
    _["metrics_numeric"] = res_numeric,
    _["metrics_categorical"] = Rcpp::List()
  );
}
