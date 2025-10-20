test_that("evaluate_imputation returns valid structure", {
  true_data <- list(
    x = c(1, 2, 3, 4, 5),
    y = c(10, 20, 30, 40, 50)
  )
  imputed_data <- list(
    x = c(1, 2, 2, 4, 6),
    y = c(10, 22, 28, 41, 49)
  )

  res <- evaluate_imputation(true_data, imputed_data, "ToyMethod")

  # --- Structural checks ---
  expect_true(is.list(res))
  expect_equal(res$method, "ToyMethod")
  expect_true(all(c("metrics_numeric", "metrics_categorical") %in% names(res)))

  # --- Numeric metrics ---
  expect_true(is.list(res$metrics_numeric))
  first_nm <- names(res$metrics_numeric[[1]])
  expect_true(all(c("RMSE", "MAE", "R2", "Correlation", "KS") %in% first_nm))

  # --- Categorical metrics (may be empty for numeric-only data) ---
  if (length(res$metrics_categorical) > 0) {
    first_cm <- names(res$metrics_categorical[[1]])
    expect_true(all(c("Accuracy", "Kappa", "MacroF1", "BalancedAccuracy") %in% first_cm))
  } else {
    expect_true(length(res$metrics_categorical) == 0)
  }
})

test_that("evaluate_imputation handles constant vectors safely", {
  true_data <- list(a = rep(5, 5))
  imputed_data <- list(a = rep(5, 5))
  res <- evaluate_imputation(true_data, imputed_data, "ConstantTest")

  expect_true(is.list(res))
  expect_equal(res$method, "ConstantTest")

  # numeric metrics should exist and RMSE = 0
  expect_true("metrics_numeric" %in% names(res))
  nm <- res$metrics_numeric[[1]]
  expect_true(is.numeric(nm$RMSE))
  expect_equal(nm$RMSE, 0)

  # Accuracy (categorical or numeric identicalness) should be 1 if present
  if (length(res$metrics_categorical) > 0) {
    cm <- res$metrics_categorical[[1]]
    expect_equal(cm$Accuracy, 1)
  }
})
