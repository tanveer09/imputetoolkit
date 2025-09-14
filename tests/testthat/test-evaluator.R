test_that("Evaluator computes metrics correctly for simple numeric data", {
  true_data <- list(
    age = c(25, 30, 40),
    income = c(50000, 60000, 70000, 8000)
  )
  imputed_data <- list(
    age = c(25, 31, 39),
    income = c(50000, 61000, 69000, 8001)
  )

  res <- evaluator(true_data, imputed_data, "mean")

  print.evaluator(res)

  # Structure
  expect_s3_class(res, "evaluator")
  expect_true("metrics" %in% names(res))
  expect_true(all(c("RMSE","MAE","R2","Correlation","KS","Accuracy") %in% names(res)))

  # Global values should be finite
  expect_true(is.finite(res$RMSE))
  expect_true(is.finite(res$MAE))
  expect_true(is.finite(res$R2))

  # Column-level metrics should exist
  expect_true("age" %in% names(res$metrics))
  expect_true("income" %in% names(res$metrics))
})

test_that("Evaluator errors on mismatched column lengths", {
  true_data <- list(
    age = c(25, 30, 40)
  )
  imputed_data <- list(
    age = c(25, 31) # length mismatch
  )

  expect_error(
    evaluator(true_data, imputed_data, "mean"),
    "same length"
  )
})
