test_that("evaluator runs with all methods on synthetic dataset", {
  # locate bundled dataset
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  expect_true(file.exists(file))

  # run evaluator (explicit filename argument)
  results <- evaluator(filename = file)

  # check that all methods are present
  expect_named(results, c("mean_mode", "median_mode", "mice"))

  # --- Mean/Mode ---
  mean_eval <- results$mean_mode
  expect_s3_class(mean_eval, "evaluator")
  expect_true(all(c("RMSE", "MAE", "R2", "Correlation", "KS", "Accuracy") %in% names(mean_eval)))

  # --- Median/Mode ---
  median_eval <- results$median_mode
  expect_s3_class(median_eval, "evaluator")
  expect_true(all(c("RMSE", "MAE", "R2", "Correlation", "KS", "Accuracy") %in% names(median_eval)))

  # --- MICE ---
  testthat::skip_if_not_installed("mice")
  mice_eval <- results$mice
  expect_s3_class(mice_eval, "evaluator")
  expect_true(all(c("RMSE", "MAE", "R2", "Correlation", "KS", "Accuracy") %in% names(mice_eval)))

  # --- Summary & print ---
  expect_s3_class(summary(mean_eval), "data.frame")
  expect_s3_class(summary(median_eval), "data.frame")

  # printing shouldn't error
  expect_invisible(print(mean_eval))
  expect_invisible(print(median_eval))

  # --- Wrapper functions ---
  metrics_df <- extract_metrics(results)
  expect_s3_class(metrics_df, "data.frame")

  # print_metrics returns a knitr_kable (visible), not invisible
  expect_s3_class(print_metrics(results), "knitr_kable")

  # plot_metrics should return a ggplot object
  expect_s3_class(plot_metrics(results, "RMSE"), "ggplot")

  # suggestion should return a character scalar
  expect_type(suggest_best_method(results, "RMSE"), "character")
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

  # Invalid metric should error
  expect_error(suggest_best_method(res, "NOT_A_METRIC"))

  # RMSE is minimization
  best_rmse <- suggest_best_method(res, "RMSE", higher_better = FALSE)
  expect_true(best_rmse %in% c("Mean/Mode", "Median/Mode", "MICE"))

  # Accuracy is maximization
  best_acc <- suggest_best_method(res, "Accuracy", higher_better = TRUE)
  expect_true(best_acc %in% c("Mean/Mode", "Median/Mode", "MICE"))
})

test_that("plot_metrics errors on invalid metric", {
  file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
  res <- evaluator(filename = file)
  expect_error(plot_metrics(res, "NotARealMetric"), "not found in data frame")
})

test_that("evaluate_imputation works with simple numeric data", {
  true_data <- list(col1 = c(1,2,3,4,5))
  imputed_data <- list(col1 = c(1,2,3,4,6))  # small error
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

