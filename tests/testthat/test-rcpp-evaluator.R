test_that("evaluate_imputation returns valid structure", {
  # create toy data
  true_data <- list(
    x = c(1, 2, 3, 4, 5),
    y = c(10, 20, 30, 40, 50)
  )
  imputed_data <- list(
    x = c(1, 2, 2, 4, 6),
    y = c(10, 22, 28, 41, 49)
  )

  res <- evaluate_imputation(true_data, imputed_data, "ToyMethod")

  #expect_s3_class(res, "list")
  expect_equal(res$method, "ToyMethod")
  expect_true("metrics" %in% names(res))
  expect_true(all(c("RMSE", "MAE", "R2", "Correlation", "KS", "Accuracy") %in% names(res)))
})
