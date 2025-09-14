#' evaluator
#'
#' Create an evaluator object for imputation quality
#'
#' @param true_data A named list of numeric vectors with actual (ground truth) values.
#' @param imputed_data A named list of numeric vectors with imputed values.
#' @param method A string giving the imputation method name (e.g., "mean", "knn").
#'
#' @return An object of class "evaluator" containing per-column and global metrics.
#' @export
#'
#' @examples
#' true_data <- list(
#'   age = c(25, 30, 40),
#'   income = c(50000, 60000, 70000)
#' )
#' imputed_data <- list(
#'   age = c(25, 31, 39),
#'   income = c(50000, 61000, 69000)
#' )
#' result <- evaluator(true_data, imputed_data, "mean")
#' print(result)
evaluator <- function(true_data, imputed_data, method) {
  res <- evaluate_imputation(true_data, imputed_data, method)
  class(res) <- "evaluator"
  res
}

#' print.evaluator
#'
#' Print method for evaluator objects
#'
#' @param x An evaluator object
#' @param ... further arguments
#' @export
print.evaluator <- function(x, ...) {
  cat("Evaluation for method:", x$method, "\n")
  cat("Global Metrics:\n")
  cat("  RMSE       :", round(x$RMSE, 4), "\n")
  cat("  MAE        :", round(x$MAE, 4), "\n")
  cat("  R^2        :", round(x$R2, 4), "\n")
  cat("  Correlation:", round(x$Correlation, 4), "\n")
  cat("  KS         :", round(x$KS, 4), "\n")
  cat("  Accuracy   :", round(x$Accuracy, 4), "\n")
  cat("\nPer-column metrics available in x$metrics\n")
}

#' summary.evaluator
#'
#' Summary method for evaluator objects
#'
#' @param object An evaluator object
#' @param ... further arguments
#' @export
summary.evaluator <- function(object, ...) {
  per_col <- do.call(rbind, lapply(object$metrics, function(m) as.data.frame(as.list(m))))
  per_col <- cbind(Column = names(object$metrics), per_col)
  rownames(per_col) <- NULL

  global <- data.frame(
    Column = "GLOBAL",
    RMSE = object$RMSE,
    MAE = object$MAE,
    R2 = object$R2,
    Correlation = object$Correlation,
    KS = object$KS,
    Accuracy = object$Accuracy
  )

  rbind(per_col, global)
}
