test_that("evaluator runs with all methods on synthetic dataset", {
  # locate bundled dataset
  file <- system.file("extdata", "synthetic_dataset.csv", package = "imputetoolkit")
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

  # evaluate_results should return invisibly
  expect_invisible(evaluate_results(results, metric = "RMSE"))
})
