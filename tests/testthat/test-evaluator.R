test_that("evaluator runs with all methods on synthetic dataset", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  expect_true(file.exists(file))

  results <- evaluator(filename = file)

  # ---- Check methods ----
  expect_named(results, c("mean_mode", "median_mode", "mice", "knn"))

  for (method_name in names(results)) {
    method_eval <- results[[method_name]]
    expect_s3_class(method_eval, "evaluator")

    # Should contain both numeric and categorical metric lists
    expect_true("metrics_numeric" %in% names(method_eval))
    expect_true("metrics_categorical" %in% names(method_eval))

    # Each numeric metric list should contain expected metrics
    if (length(method_eval$metrics_numeric) > 0) {
      nm <- names(method_eval$metrics_numeric[[1]])
      expect_true(all(c("RMSE", "MAE", "R2", "Correlation", "KS") %in% nm))
    }

    # Each categorical metric list should contain expected metrics
    if (length(method_eval$metrics_categorical) > 0) {
      cm <- names(method_eval$metrics_categorical[[1]])
      expect_true(all(c("Accuracy", "Kappa", "MacroF1", "BalancedAccuracy") %in% cm))
    }
  }

  # ---- Summary & print ----
  s <- summary(results$mean_mode)
  expect_true(is.list(s))
  expect_true(all(c("numeric_metrics", "categorical_metrics") %in% names(s)))
  expect_s3_class(print(results$mean_mode), "evaluator")

  # ---- Wrapper functions ----
  metrics_df <- extract_metrics(results)
  expect_s3_class(metrics_df, "data.frame")
  expect_true(any(grepl("RMSE|MAE|Accuracy", colnames(metrics_df))))

  expect_s3_class(print_metrics(results), "knitr_kable")
  expect_s3_class(plot_metrics(results, "RMSE"), "ggplot")

  # ---- Suggest best ----
  best_rmse <- suggest_best_method(results, "RMSE")
  expect_true(is.list(best_rmse))
  expect_true(best_rmse$numeric %in% c("Mean/Mode", "Median/Mode", "MICE", "KNN"))

  best_all <- suggest_best_method(results, "ALL")
  expect_true(is.list(best_all))
  all_metrics <- unlist(c(best_all$numeric, best_all$categorical), use.names = FALSE)
    expect_true(all(all_metrics %in% c(
      "RMSE", "MAE", "R2", "Correlation", "KS",
      "Accuracy", "Kappa", "F1", "MacroF1", "BalancedAccuracy")))

})

test_that("evaluator fails with invalid inputs", {
  expect_error(evaluator(), "Please provide either a filename or a data.frame")
  expect_error(evaluator(filename = "nonexistent.csv"))

  tmp <- tempfile(fileext = ".xyz")
  writeLines("dummy", tmp)
  expect_error(evaluator(filename = tmp), "Unsupported file type")
})

test_that("evaluator is reproducible with fixed seed", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  res1 <- evaluator(filename = file)
  res2 <- evaluator(filename = file)
  expect_equal(res1$mean_mode$metrics_numeric[[1]], res2$mean_mode$metrics_numeric[[1]])
})

test_that("extract_metrics validates inputs", {
  expect_error(extract_metrics("not a list"), "res must be a list")
})

test_that("suggest_best_method works correctly", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  res <- evaluator(filename = file)

  expect_error(suggest_best_method(res, "NOT_A_METRIC"))

  for (m in c("RMSE", "MAE", "R2", "KS", "Accuracy")) {
    best <- suggest_best_method(res, m)
    expect_true(is.list(best))
    if (m %in% c("RMSE", "MAE", "R2", "KS")) {
      expect_true(best$numeric %in% c("Mean/Mode", "Median/Mode", "MICE", "KNN"))
    } else {
      expect_true(best$categorical %in% c("Mean/Mode", "Median/Mode", "MICE", "KNN"))
    }
  }
})

test_that("plot_metrics works correctly", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  res <- evaluator(filename = file)

  expect_error(plot_metrics(res, "NotARealMetric"), "not found")

  expect_s3_class(plot_metrics(res, "RMSE"), "ggplot")
  expect_s3_class(plot_metrics(res, "ALL"), "ggplot")
})

test_that("evaluate_imputation works with simple numeric data", {
  true_data    <- list(col1 = c(1,2,3,4,5))
  imputed_data <- list(col1 = c(1,2,3,4,6))

  res <- evaluate_imputation(true_data, imputed_data, "toy_method")

  expect_true(is.list(res))
  expect_equal(res$method, "toy_method")
  expect_true(all(c("metrics_numeric","metrics_categorical") %in% names(res)))

  # numeric metrics exist & look reasonable
  expect_true(length(res$metrics_numeric) >= 1)
  nm <- res$metrics_numeric[[1]]
  expect_true(all(c("RMSE","MAE","R2","Correlation","KS","Accuracy") %in% names(nm)))
  expect_gte(nm$RMSE, 0)

  # categorical should be empty here
  expect_equal(length(res$metrics_categorical), 0)
})

test_that("print and summary invisibility", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  res <- evaluator(filename = file)
  expect_invisible(print(res$mean_mode))
  s <- summary(res$mean_mode)
  expect_true(is.list(s))
  expect_true(all(c("numeric_metrics","categorical_metrics") %in% names(s)))
})

# --- Additional edge tests for KNN ---
test_that("KNN imputation produces reasonable values", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  res <- evaluator(filename = file)
  knn_eval <- res$knn
  expect_s3_class(knn_eval, "evaluator")

  if (length(knn_eval$metrics_numeric) > 0) {
    vals <- unlist(lapply(knn_eval$metrics_numeric, function(x) x$RMSE))
    expect_true(all(is.numeric(vals)))
    expect_true(all(vals >= 0))
  }
})
