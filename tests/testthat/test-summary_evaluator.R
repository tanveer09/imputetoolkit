test_that("Summary method returns combined data.frame", {
  true_data <- list(
    x = c(1,2,3),
    y = c(10,20,30, 32, 33, 35, 35)
  )
  imputed_data <- list(
    x = c(1,2,3),
    y = c(10,20,30, 32, 33, 35, 35)
  )

  res <- evaluator(true_data, imputed_data, "knn")
  sm <- summary(res)

  print.evaluator(res)

  expect_s3_class(sm, "data.frame")
  expect_true("GLOBAL" %in% sm$Column)
  expect_true(all(c("RMSE","MAE","R2","Correlation","KS","Accuracy") %in% names(sm)))
})

