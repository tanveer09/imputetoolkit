#' evaluator
#'
#' Full pipeline: data cleaning, missingness injection, imputation, and evaluation.
#'
#' @param data A data.frame that has already been loaded in R. Optional if `filename` is provided.
#' @param filename Path to a data file (CSV, TSV, TXT, XLSX, XLS, RDS). Optional if `data` is provided.

#' @return A list of evaluator objects (one per imputation method).
#' @export
#'
#' @importFrom utils read.csv read.delim
#' @importFrom readxl read_excel
#' @importFrom stats median
#' @importFrom mice mice complete
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr all_of

#'
#' @examples
#' # Example 1: Using a file shipped with the package
#' file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
#' \dontrun{
#' res <- evaluator(filename = file)
#' print(res$mean_mode)
#' summary(res$median_mode)
#' }
#'
#' # Example 2: Passing a pre-loaded data.frame
#' df <- utils::read.csv(system.file("extdata", "sample_dataset.csv",
#'                                   package = "imputetoolkit"),
#'                       stringsAsFactors = TRUE)
#' \dontrun{
#' res <- evaluator(data = df)
#' print(res$mean_mode)
#' }
evaluator <- function(data = NULL, filename = NULL) {
  set.seed(123) # reproducibility

  # ---- 0. Load data ----
  if (!is.null(data)) {
    raw_data <- data
  } else if (!is.null(filename)) {
    file_ext <- tools::file_ext(filename)

    if (file_ext == "csv") {
      raw_data <- utils::read.csv(filename, stringsAsFactors = TRUE)
    } else if (file_ext %in% c("tsv", "txt")) {
      raw_data <- utils::read.delim(filename, stringsAsFactors = TRUE)
    } else if (file_ext %in% c("xlsx", "xls")) {
      raw_data <- readxl::read_excel(filename)
    } else if (file_ext %in% c("rds")) {
      raw_data <- readRDS(filename)
    } else {
      stop("Unsupported file type: ", file_ext)
    }
  } else {
    stop("Please provide either a filename or a data.frame.")
  }

  # Convert empty strings to NA for categorical cols
  for (col in names(raw_data)) {
    if (is.factor(raw_data[[col]]) || is.character(raw_data[[col]])) {
      raw_data[[col]][raw_data[[col]] == ""] <- NA
    }
  }

  # ---- 1. Inject extra missingness ----
  raw_data_col_missingness <- colSums(is.na(raw_data)) / nrow(raw_data) * 100
  raw_data_col_missingness <- raw_data_col_missingness[raw_data_col_missingness > 0]

  raw_data_selected_indices <- lapply(names(raw_data_col_missingness), function(col) {
    non_missing_idx <- which(!is.na(raw_data[[col]]))
    perc <- raw_data_col_missingness[col]
    n_to_sample <- ceiling(length(non_missing_idx) * perc / 100)
    sample(non_missing_idx, n_to_sample)
  })
  names(raw_data_selected_indices) <- names(raw_data_col_missingness)

  raw_data_modified <- raw_data
  for (col in names(raw_data_selected_indices)) {
    idx <- raw_data_selected_indices[[col]]
    raw_data_modified[idx, col] <- NA
  }
  # Re-enforce factors after NA injection
  cat_cols <- names(raw_data)[sapply(raw_data, is.factor)]
  raw_data_modified[cat_cols] <- lapply(raw_data_modified[cat_cols], factor)

  # ---- 2. Get mode of the given column ----
  get_mode <- function(x) {
    ux <- unique(x[!is.na(x)])
    ux[which.max(tabulate(match(x, ux)))]
  }

  # ---- 3. Imputation methods ----
  # Mean/Mode
  data_mean <- raw_data_modified
  for (col in names(data_mean)) {
    if (is.numeric(data_mean[[col]])) {
      val <- mean(raw_data[[col]], na.rm = TRUE)
      data_mean[[col]][is.na(data_mean[[col]])] <- val
    } else {
      val <- get_mode(raw_data[[col]])
      data_mean[[col]][is.na(data_mean[[col]])] <- val
    }
  }

  # Median/Mode
  data_median <- raw_data_modified
  for (col in names(data_median)) {
    if (is.numeric(data_median[[col]])) {
      val <- stats::median(raw_data[[col]], na.rm = TRUE)
      data_median[[col]][is.na(data_median[[col]])] <- val
    } else {
      val <- get_mode(raw_data[[col]])
      data_median[[col]][is.na(data_median[[col]])] <- val
    }
  }

  # MICE (Predictive Mean Matching)
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' needed for this function to work. Please install it.")
  }
  mice_imp <- mice::mice(raw_data_modified, m = 1, method = "pmm", maxit = 5, seed = 123)
  data_mice <- mice::complete(mice_imp)

  # ---- 4. Collect true vs. imputed values for evaluation ----
  build_eval_data <- function(imputed_data) {
    true_list <- list()
    imp_list  <- list()
    for (col in names(raw_data_selected_indices)) {
      idx <- raw_data_selected_indices[[col]]
      true_list[[col]] <- raw_data[[col]][idx]
      imp_list[[col]]  <- imputed_data[[col]][idx]
    }
    list(true = true_list, imputed = imp_list)
  }

  eval_mean   <- build_eval_data(data_mean)
  eval_median <- build_eval_data(data_median)
  eval_mice   <- build_eval_data(data_mice)

  # ---- 5. Call Rcpp evaluator ----
  res_mean   <- evaluate_imputation(eval_mean$true,   eval_mean$imputed,   "Mean/Mode")
  res_median <- evaluate_imputation(eval_median$true, eval_median$imputed, "Median/Mode")
  res_mice   <- evaluate_imputation(eval_mice$true,   eval_mice$imputed,   "MICE")

  class(res_mean)   <- "evaluator"
  class(res_median) <- "evaluator"
  class(res_mice)   <- "evaluator"

  # Return all methods
  list(
    mean_mode   = res_mean,
    median_mode = res_median,
    mice        = res_mice
  )
}


#' print.evaluator
#'
#' Print method for evaluator objects
#'
#' @param x An evaluator object
#' @param ... Additional arguments (ignored)
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
#' @param ... Additional arguments (ignored)
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



#' extract_metrics
#'
#' Extract evaluation metrics from Rcpp evaluator output
#'
#' @param res The result object returned by evaluator() (a list of methods)
#' @return A data frame with evaluation metrics for each method
#' @keywords internal
#' @noRd
extract_metrics <- function(res) {
  if (!is.list(res)) stop("res must be a list returned by evaluator().")

  metrics_list <- lapply(names(res), function(method_name) {
    method_res <- res[[method_name]]

    data.frame(
      Method      = method_res$method,
      RMSE        = method_res$RMSE,
      MAE         = method_res$MAE,
      R2          = method_res$R2,
      Correlation = method_res$Correlation,
      KS          = method_res$KS,
      Accuracy    = method_res$Accuracy,
      stringsAsFactors = FALSE
    )
  })

  # Base R only, no dplyr
  metrics_df <- do.call(rbind, metrics_list)
  rownames(metrics_df) <- NULL
  return(metrics_df)
}



#' Print evaluation metrics for imputation methods
#'
#' @param x A result list returned by evaluator()
#' @return Prints the metrics in a clean table format
#' @export
print_metrics <- function(x) {

  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required. Please install it.")
  }
  metrics_df <- extract_metrics(x)

  knitr::kable(
    metrics_df,
    digits = 4,
    caption = "Comparison of Imputation Methods"
  )
}


#' Plot evaluation metrics
#'
#' @param x A result list returned by evaluator()
#' @param metric Character, which metric to plot (e.g., "RMSE", "MAE", "R2", "KS", "Accuracy", or "ALL")
#' @return A ggplot object
#' @export
plot_metrics <- function(x, metric = "RMSE") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Please install it.")
  }

  metrics_df <- extract_metrics(x)
  valid_metrics <- c("RMSE", "MAE", "R2", "KS", "Accuracy")

  if (metric == "ALL") {
    long_df <- tidyr::pivot_longer(
      metrics_df,
      cols = dplyr::all_of(valid_metrics),
      names_to = "Metric",
      values_to = "Value"
    )

    p <- ggplot2::ggplot(long_df, ggplot2::aes(x = .data$Method, y = .data$Value, fill = .data$Method)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::facet_wrap(~ Metric, scales = "free_y") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Comparison of Methods across Metrics",
        y = "Value", x = "Method"
      ) +
      ggplot2::theme(legend.position = "none")

    return(p)

  } else {
    if (!(metric %in% valid_metrics)) {
      stop("Metric ", metric, " not found. Please use one of: ", paste(valid_metrics, collapse = ", "), " or ALL")
    }

    p <- ggplot2::ggplot(metrics_df, ggplot2::aes(x = .data$Method, y = .data[[metric]], fill = .data$Method)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Comparison of Methods on", metric),
        y = metric, x = "Method"
      ) +
      ggplot2::theme(legend.position = "none")

    return(p)
  }
}


#' Suggest the best imputation method
#'
#' @param x A result list returned by evaluator()
#' @param metric Character, metric to optimize
#'   ("RMSE", "MAE", "R2", "KS", "Accuracy", or "ALL").
#'   If "ALL" or missing, the function compares across all metrics.
#' @return If a single metric is provided, returns the best method (character).
#'   If "ALL", returns a named list showing the metrics associated with each best method.
#' @export
suggest_best_method <- function(x, metric = "ALL") {

  metrics_df <- extract_metrics(x)

  higher_metrics <- c("R2", "KS", "Accuracy")
  lower_metrics  <- c("RMSE", "MAE")
  valid_metrics  <- c(higher_metrics, lower_metrics)

  if (metric == "ALL" || is.null(metric)) {
    best_methods <- list()

    for (m in valid_metrics) {
      if (!(m %in% colnames(metrics_df))) next

      if (m %in% higher_metrics) {
        best_idx <- which.max(metrics_df[[m]])
      } else {
        best_idx <- which.min(metrics_df[[m]])
      }
      method <- metrics_df$Method[best_idx]
      best_methods[[m]] <- method
    }

    # Group metrics by best method
    grouped <- split(names(best_methods), unlist(best_methods))

    # Build summary message
    msg <- paste(
      vapply(names(grouped),
             function(method) {
               metrics <- paste(grouped[[method]], collapse = ", ")
               paste0("    As per ", metrics, " metrics: ", method)
             },
             character(1L)),
      collapse = "\n"
    )

    message("Suggested best imputation methods across metrics:\n", msg)
    return(grouped)

  } else {
    if (!(metric %in% valid_metrics)) {
      stop("Unknown metric: ", metric,
           ". Please use one of: ", paste(valid_metrics, collapse = ", "), ", or ALL")
    }

    if (metric %in% higher_metrics) {
      best_idx <- which.max(metrics_df[[metric]])
    } else {
      best_idx <- which.min(metrics_df[[metric]])
    }

    best_method <- metrics_df$Method[best_idx]
    message("Suggested best imputation method based on ", metric, ": ", best_method)
    return(best_method)
  }
}



#' Evaluate results: print, plot, and suggest best method
#'
#' @param res The result object returned from evaluator()
#' @param metric Character, which metric to optimize (default = "RMSE")
#' @return Prints a table, plots results, and returns the suggested method
#' @export
evaluate_results <- function(res, metric = "RMSE") {
  metrics_df <- extract_metrics(res)
  print_metrics(metrics_df)
  print(plot_metrics(metrics_df, metric))
  best <- suggest_best_method(metrics_df, metric)
  invisible(best)
}

