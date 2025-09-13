#' Create an evaluator object
#'
#' @param true_vals Numeric vector of actual values (ground truth).
#' @param imputed_vals Numeric vector of imputed values (same length as true_vals).
#' @param method Character string naming the imputation method.
#'
#' @return An object of class "evaluator" containing method name and metrics.
#' @export
#'
#' @examples
#' true_vals <- c(1, 2, 3, 4, 5)
#' imputed_vals <- c(1.1, 1.9, 3.2, 3.8, 5.1)
#' result <- evaluator(true_vals, imputed_vals, "mean")
#' print(result)
evaluator <- function(true_vals, imputed_vals, method) {
  res <- evaluate_imputation(true_vals, imputed_vals, method)
  class(res) <- "evaluator"
  res
}

#' Print method for evaluator objects
#' @param x An evaluator object
#' @param ... further arguments
#' @export
print.evaluator <- function(x, ...) {
  cat("Evaluation for method:", x$method, "\n")
  cat("  RMSE:", round(x$RMSE, 4), "\n")
  cat("  MAE :", round(x$MAE, 4), "\n")
  cat("  R^2 :", round(x$R2, 4), "\n")
}

#' Summary method for evaluator objects
#' @param object An evaluator object
#' @param ... further arguments
#' @export
summary.evaluator <- function(object, ...) {
  data.frame(
    Method = object$method,
    RMSE   = object$RMSE,
    MAE    = object$MAE,
    R2     = object$R2
  )
}
