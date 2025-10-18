test_that("evaluator runs with all methods on synthetic dataset", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  expect_true(file.exists(file))

  results <- evaluator(filename = file)

  # ---- Check methods ----
  expect_named(results, c("mean_mode", "median_mode", "mice", "knn"))

  for (method_name in names(results)) {
    method_eval <- results[[method_name]]
    expect_s3_class(method_eval, "evaluator")
    expect_true(all(c("RMSE", "MAE", "R2", "Correlation", "KS", "Accuracy") %in% names(method_eval)))
  }

  # ---- Summary & print ----
  expect_s3_class(summary(results$mean_mode), "data.frame")
  expect_invisible(print(results$mean_mode))

  # ---- Wrapper functions ----
  metrics_df <- extract_metrics(results)
  expect_s3_class(metrics_df, "data.frame")
  expect_true(all(c("RMSE", "MAE", "R2", "Correlation", "KS", "Accuracy") %in% colnames(metrics_df)))

  expect_s3_class(print_metrics(results), "knitr_kable")
  expect_s3_class(plot_metrics(results, "RMSE"), "ggplot")

  # ---- Suggest best ----
  best_rmse <- suggest_best_method(results, "RMSE")
  expect_type(best_rmse, "character")
  expect_true(best_rmse %in% c("Mean/Mode", "Median/Mode", "MICE", "KNN"))

  best_all <- suggest_best_method(results, "ALL")
  expect_type(best_all, "list")
  expect_true(all(names(best_all) %in% c("Mean/Mode", "Median/Mode", "MICE", "KNN")))
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
  expect_equal(res1$mean_mode$RMSE, res2$mean_mode$RMSE)
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
    expect_true(best %in% c("Mean/Mode", "Median/Mode", "MICE", "KNN"))
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
  true_data <- list(col1 = c(1,2,3,4,5))
  imputed_data <- list(col1 = c(1,2,3,4,6))
  res <- evaluate_imputation(true_data, imputed_data, "toy_method")

  expect_true(is.list(res))
  expect_true(all(c("RMSE","MAE","R2","Correlation","KS","Accuracy") %in% names(res)))
  expect_equal(res$method, "toy_method")
})

test_that("print and summary invisibility", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  res <- evaluator(filename = file)
  expect_invisible(print(res$mean_mode))
  expect_s3_class(summary(res$mean_mode), "data.frame")
})

# --- Additional edge tests for KNN ---
test_that("KNN imputation produces reasonable values", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  res <- evaluator(filename = file)
  knn_eval <- res$knn
  expect_s3_class(knn_eval, "evaluator")
  expect_true(is.numeric(knn_eval$RMSE))
  expect_true(knn_eval$RMSE >= 0)
})
