test_that("print.evaluator outputs expected text", {
  true_data <- list(
    x = c(1, 2, 3),
    y = c(10, 20, 30)
  )
  imputed_data <- list(
    x = c(1, 2, 4),
    y = c(10, 21, 29)
  )

  res <- evaluator(true_data, imputed_data, "mean")

  # Capture output of print()
  expect_output(print(res), "Evaluation for method: mean")
  expect_output(print(res), "Global Metrics:")
  expect_output(print(res), "RMSE")
  expect_output(print(res), "MAE")
  expect_output(print(res), "R\\^2")   # regex escape for caret
  expect_output(print(res), "Correlation")
  expect_output(print(res), "KS")
  expect_output(print(res), "Accuracy")
  expect_output(print(res), "Per-column metrics available")
})

