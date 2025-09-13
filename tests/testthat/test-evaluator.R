test_that("Evaluator computes metrics correctly for mean imputation", {
  true_vals <- c(1, 2, 3, 4, 5)
  imputed_vals <- c(1.1, 1.9, 3.2, 3.8, 5.1)

  res <- evaluator(true_vals, imputed_vals, "mean")

  expect_s3_class(res, "evaluator")
  expect_named(res, c("method", "RMSE", "MAE", "R2"))

  # Approximate expectations
  expect_equal(res$method, "mean")
  expect_equal(round(res$RMSE, 3), 0.141, tolerance = 1e-3)
  expect_equal(round(res$MAE, 3), 0.120, tolerance = 1e-3)
  expect_true(res$R2 > 0.95)
})

test_that("Evaluator throws error when input lengths mismatch", {
  true_vals <- c(1, 2, 3, 4, 5)
  imputed_vals <- c(1, 2, 3)  # shorter

  expect_error(
    evaluator(true_vals, imputed_vals, "mean"),
    "must have the same length"
  )
})
