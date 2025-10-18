#' evaluator
#'
#' Full pipeline: data cleaning, missingness injection, imputation, and evaluation.
#'
#' @param data A data.frame that has already been loaded in R. Optional if `filename` is provided.
#' @param filename Path to a data file (CSV, TSV, TXT, XLSX, XLS, RDS). Optional if `data` is provided.

#' @return A list of evaluator objects (one per imputation method).
#' @seealso vignette("imputetoolkit")
#' For a complete tutorial, see the package vignette:
#' \code{vignette("imputetoolkit")}
#' @export
#'
#' @importFrom utils read.csv read.delim
#' @importFrom readxl read_excel
#' @importFrom stats median
#' @importFrom mice mice complete
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr all_of
#' @importFrom foreach foreach "%dopar%"
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom stats complete.cases

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

  ##########################
  # ---- 0. Load data ---- #
  ##########################

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
  #########################################
  # ---- 1. Inject extra missingness ---- #
  #########################################
  # --- Identify columns with missing values ---
  raw_data_col_missingness <- colSums(is.na(raw_data)) / nrow(raw_data) * 100
  raw_data_col_missingness <- raw_data_col_missingness[raw_data_col_missingness > 0]

  # --- Filter numeric columns among those with missing values ---
  numeric_missing_cols <- names(raw_data_col_missingness)[
    sapply(raw_data[names(raw_data_col_missingness)], is.numeric)
  ]

  # --- Compute min and max for each numeric column with missing values ---
  min_max_missing <- lapply(numeric_missing_cols, function(col) {
    vals <- raw_data[[col]]
    c(min = min(vals, na.rm = TRUE), max = max(vals, na.rm = TRUE))
  })

  # --- Convert to a data frame for easy reference ---
  min_max_missing <- as.data.frame(do.call(rbind, min_max_missing))
  rownames(min_max_missing) <- numeric_missing_cols

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

  #############################################
  # ---- 2. Get mode of the given column ---- #
  #############################################
  get_mode <- function(x) {
    ux <- unique(x[!is.na(x)])
    ux[which.max(tabulate(match(x, ux)))]
  }

  ###################################
  # ---- 3. Imputation methods ---- #
  ###################################

  # ---- Mean/Mode Imputation with Skewness Check ----
  data_mean <- raw_data_modified

  for (col in names(data_mean)) {
    if (is.numeric(data_mean[[col]])) {

      # 1 Calculate skewness (ignoring missing values)
      skew_val <- e1071::skewness(raw_data[[col]], na.rm = TRUE)

      # 2 If highly skewed (|skewness| > 1), apply log transform before computing mean
      if (abs(skew_val) > 1) {
        # Add a small constant to handle zeros or negatives
        temp <- raw_data[[col]]
        min_val <- min(temp, na.rm = TRUE)
        if (min_val <= 0) {
          temp <- temp + abs(min_val) + 1
        }

        # Apply log transformation
        log_vals <- log(temp)

        # Compute mean in log-space, then back-transform
        val <- exp(mean(log_vals, na.rm = TRUE))  # geometric mean
      } else {
        # 3 If not skewed, just use the arithmetic mean
        val <- mean(raw_data[[col]], na.rm = TRUE)
      }

      # 4 Replace missing values with computed value
      data_mean[[col]][is.na(data_mean[[col]])] <- val

    } else {
      # 5 For categorical columns, use mode
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


  ##############################################
  # ---- KNN Imputation ---- #
  ##############################################
  message("Running mixed-type parallel KNN imputation using FNN...")

  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' required for fast KNN imputation. Please install it.")
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' required for parallel execution. Please install it.")
  }
  if (!requireNamespace("foreach", quietly = TRUE)) stop("Package 'foreach' required.")

  data_knn <- raw_data_modified

  # Separate numeric and categorical columns
  num_cols <- names(which(sapply(data_knn, is.numeric)))
  cat_cols <- names(which(!sapply(data_knn, is.numeric)))

  # Standardize numeric data
  num_data <- data_knn[, num_cols, drop = FALSE]
  num_data_scaled <- scale(num_data)

  # Prepare predictor matrix (numeric only)
  pred_data <- as.data.frame(num_data_scaled)

  # ---- Setup parallel backend safely ----
  total_cores <- parallel::detectCores(logical = TRUE)

  # detect if under R CMD check / testthat
  is_check_env <- isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE"))) ||
    ("CHECK" %in% toupper(names(Sys.getenv()))) ||
    ("TESTTHAT" %in% toupper(names(Sys.getenv())))

  # Use 1 core inside check/test environments to avoid socket issues
  cores <- if (is_check_env) 1 else max(1, total_cores - 1)

  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    message("Using ", cores, " cores for KNN imputation...")
  } else {
    foreach::registerDoSEQ()
    message("Using sequential mode for KNN imputation (testing environment)...")
  }

  # ---- Impute numeric columns ----
  imputed_numeric <- foreach::foreach(col = num_cols, .packages = "FNN") %dopar% {
    col_data <- num_data_scaled[, col]
    na_idx <- which(is.na(col_data))
    if (length(na_idx) == 0) return(col_data)

    pred_cols <- setdiff(num_cols, col)
    if (length(pred_cols) == 0) return(col_data)

    train_idx <- which(!is.na(col_data) & complete.cases(pred_data[, pred_cols, drop = FALSE]))
    test_idx  <- which(is.na(col_data) & complete.cases(pred_data[, pred_cols, drop = FALSE]))

    if (length(train_idx) == 0 || length(test_idx) == 0) return(col_data)

    train_x <- pred_data[train_idx, pred_cols, drop = FALSE]
    test_x  <- pred_data[test_idx, pred_cols, drop = FALSE]

    # Find nearest neighbors
    knn_res <- FNN::get.knnx(train_x, test_x, k = 3)

    # Impute by mean of neighbor values
    imp_vals <- vapply(seq_len(nrow(knn_res$nn.index)), function(i) {
      neighbor_idx <- knn_res$nn.index[i, ]
      mean(col_data[train_idx[neighbor_idx]], na.rm = TRUE)
    }, numeric(1))

    col_data[test_idx] <- imp_vals
    col_data
  }

  # ---- Impute categorical columns ----
  imputed_categorical <- foreach::foreach(col = cat_cols, .packages = "FNN") %dopar% {
    col_data <- data_knn[[col]]
    na_idx <- which(is.na(col_data))
    if (length(na_idx) == 0) return(as.character(col_data))  # Always return character

    # Use numeric predictors for neighbor search
    if (length(num_cols) == 0) return(as.character(col_data))

    train_idx <- which(!is.na(col_data) & complete.cases(pred_data))
    test_idx  <- which(is.na(col_data) & complete.cases(pred_data))
    if (length(train_idx) == 0 || length(test_idx) == 0)
      return(as.character(col_data))

    train_x <- pred_data[train_idx, , drop = FALSE]
    test_x  <- pred_data[test_idx, , drop = FALSE]

    # Nearest numeric neighbors
    knn_res <- FNN::get.knnx(train_x, test_x, k = 3)

    # Safe mode finder
    get_mode <- function(x) {
      ux <- unique(x[!is.na(x)])
      if (length(ux) == 0) return(NA_character_)
      ux[which.max(tabulate(match(x, ux)))]
    }

    imp_vals <- vapply(seq_len(nrow(knn_res$nn.index)), function(i) {
      neighbor_idx <- knn_res$nn.index[i, ]
      get_mode(as.character(col_data[train_idx[neighbor_idx]]))
    }, character(1))

    # Replace missing with imputed
    col_data[test_idx] <- imp_vals
    as.character(col_data)  # Always return character
  }

  # Safely stop the cluster only if it exists
  if (exists("cl") && inherits(cl, "cluster")) {
    parallel::stopCluster(cl)
    message("Parallel mixed-type KNN imputation complete (parallel).")
  } else {
    message("KNN imputation completed in sequential mode (no cluster to stop).")
  }

  # ---- Merge results ----
  for (i in seq_along(num_cols)) {
    num_data_scaled[, num_cols[i]] <- imputed_numeric[[i]]
  }
  num_data_imputed <- sweep(num_data_scaled, 2, attr(num_data_scaled, "scaled:scale"), `*`)
  num_data_imputed <- sweep(num_data_imputed, 2, attr(num_data_scaled, "scaled:center"), `+`)
  data_knn[num_cols] <- num_data_imputed

  # Convert back categorical to proper factor
  for (i in seq_along(cat_cols)) {
    data_knn[[cat_cols[i]]] <- factor(imputed_categorical[[i]], levels = levels(raw_data[[cat_cols[i]]]))
  }


  # ---- Final fallback for any remaining NAs ----
  for (col in names(data_knn)) {
    if (anyNA(data_knn[[col]])) {
      if (is.numeric(data_knn[[col]])) {
        data_knn[[col]][is.na(data_knn[[col]])] <- mean(data_knn[[col]], na.rm = TRUE)
      } else {
        ux <- unique(data_knn[[col]][!is.na(data_knn[[col]])])
        if (length(ux) > 0) {
          mode_val <- ux[which.max(tabulate(match(data_knn[[col]], ux)))]
          data_knn[[col]][is.na(data_knn[[col]])] <- mode_val
        }
      }
    }
  }


  ###############################################################
  # ---- 4. Collect true vs. imputed values for evaluation ---- #
  ###############################################################

  # --- Collects true & imputed, then scales them ---
  build_eval_data <- function(imputed_data) {
    true_list <- list()
    imp_list  <- list()

    for (col in names(raw_data_selected_indices)) {
      idx <- raw_data_selected_indices[[col]]
      true_vals <- raw_data[[col]][idx]
      imp_vals  <- imputed_data[[col]][idx]

      # Apply scaling only if column has recorded min-max
      if (col %in% rownames(min_max_missing)) {
        min_val <- min_max_missing[col, "min"]
        max_val <- min_max_missing[col, "max"]

        # Avoid division by zero if all values are the same
        if (max_val > min_val) {
          true_vals <- (true_vals - min_val) / (max_val - min_val)
          imp_vals  <- (imp_vals  - min_val) / (max_val - min_val)
        }
      }

      true_list[[col]] <- true_vals
      imp_list[[col]]  <- imp_vals
    }

    list(true = true_list, imputed = imp_list)
  }


  eval_mean   <- build_eval_data(data_mean)
  eval_median <- build_eval_data(data_median)
  eval_mice   <- build_eval_data(data_mice)
  eval_knn <- build_eval_data(data_knn)


  # ---- 5. Call Rcpp evaluator ----
  res_mean   <- evaluate_imputation(eval_mean$true,   eval_mean$imputed,   "Mean/Mode")
  res_median <- evaluate_imputation(eval_median$true, eval_median$imputed, "Median/Mode")
  res_mice   <- evaluate_imputation(eval_mice$true,   eval_mice$imputed,   "MICE")
  res_knn <- evaluate_imputation(eval_knn$true, eval_knn$imputed, "KNN")

  class(res_mean)   <- "evaluator"
  class(res_median) <- "evaluator"
  class(res_mice)   <- "evaluator"
  class(res_knn) <- "evaluator"

  # Attach the eval pairs so plotting functions can access them later
  res_mean$eval   <- eval_mean
  res_median$eval <- eval_median
  res_mice$eval   <- eval_mice
  res_knn$eval    <- eval_knn

  # Return all methods
  list(
    mean_mode   = res_mean,
    median_mode = res_median,
    mice        = res_mice,
    knn         = res_knn
  )

}


#' print.evaluator
#'
#' Print method for evaluator objects
#'
#' @param x An evaluator object
#' @param ... Additional arguments (ignored)
#' @seealso vignette("imputetoolkit")
#' For a complete tutorial, see the package vignette:
#' \code{vignette("imputetoolkit")}
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
#' @seealso vignette("imputetoolkit")
#' For a complete tutorial, see the package vignette:
#' \code{vignette("imputetoolkit")}
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
#' @seealso vignette("imputetoolkit")
#' For a complete tutorial, see the package vignette:
#' \code{vignette("imputetoolkit")}
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
#' @seealso vignette("imputetoolkit")
#' For a complete tutorial, see the package vignette:
#' \code{vignette("imputetoolkit")}
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


#' Plot density comparisons of true vs imputed values for one column
#'
#' @param eval_list A named list containing evaluation data for each method
#'   (each should have $true and $imputed lists)
#' @param col_name Character, name of the column to plot
#' @return A ggplot2 object showing density overlays of true vs imputed data
#' @export
plot_density_per_column <- function(eval_list, col_name) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }

  # Prepare data for plotting
  plot_data <- lapply(names(eval_list), function(method) {
    eval_data <- eval_list[[method]]
    if (is.null(eval_data$true[[col_name]])) return(NULL)

    data.frame(
      Value = c(eval_data$true[[col_name]], eval_data$imputed[[col_name]]),
      Type = rep(c("True", "Imputed"), each = length(eval_data$true[[col_name]])),
      Method = method
    )
  })

  plot_data <- do.call(rbind, plot_data)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, fill = Type)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::facet_wrap(~ Method, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Density Comparison for", col_name),
      x = col_name, y = "Density"
    ) +
    ggplot2::theme(legend.position = "top")
}


#' Plot combined density comparisons for all columns and methods
#'
#' @param eval_list A named list containing evaluation data for each method
#'   (mean_mode, median_mode, mice, knn)
#' @return A faceted ggplot object showing all density comparisons
#' @export
plot_density_all <- function(eval_list) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }

  # Combine all true vs imputed data for every column and method
  all_data <- lapply(names(eval_list), function(method) {
    eval_data <- eval_list[[method]]
    do.call(rbind, lapply(names(eval_data$true), function(col) {
      data.frame(
        Value = c(eval_data$true[[col]], eval_data$imputed[[col]]),
        Type = rep(c("True", "Imputed"), each = length(eval_data$true[[col]])),
        Method = method,
        Column = col
      )
    }))
  })

  all_data <- do.call(rbind, all_data)

  ggplot2::ggplot(all_data, ggplot2::aes(x = Value, fill = Type)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::facet_grid(Column ~ Method, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Density Comparison Across All Columns and Methods",
      x = "Value", y = "Density"
    ) +
    ggplot2::theme(legend.position = "top")
}

#' Suggest the best imputation method
#'
#' @param x A result list returned by evaluator()
#' @param metric Character, metric to optimize
#'   ("RMSE", "MAE", "R2", "KS", "Accuracy", or "ALL").
#'   If "ALL" or missing, the function compares across all metrics.
#' @return If a single metric is provided, returns the best method (character).
#'   If "ALL", returns a named list showing the metrics associated with each best method.
#' @seealso vignette("imputetoolkit")
#' For a complete tutorial, see the package vignette:
#' \code{vignette("imputetoolkit")}
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
#' @seealso vignette("imputetoolkit")
#' For a complete tutorial, see the package vignette:
#' \code{vignette("imputetoolkit")}
#' @export
evaluate_results <- function(res, metric = "RMSE") {
  print_metrics(res)
  plot_metrics(res, metric)
  best <- suggest_best_method(res, metric)
  invisible(best)
}

#' Collect true-vs-imputed evaluation lists from evaluator() output
#' @param res The list returned by evaluator()
#' @return A named list: each element has $true and $imputed lists by column
#' @export
get_eval_list <- function(res) {
  stopifnot(is.list(res))
  out <- list()
  for (m in names(res)) {
    if (!is.null(res[[m]]$eval) &&
        is.list(res[[m]]$eval) &&
        !is.null(res[[m]]$eval$true) &&
        !is.null(res[[m]]$eval$imputed)) {
      out[[m]] <- res[[m]]$eval
    }
  }
  if (length(out) == 0) {
    stop("No embedded eval data found. Update evaluator() to attach $eval per method.")
  }
  out
}

#' Plot density comparisons of true vs imputed values for one column
#' @param eval_list Output of get_eval_list(res)
#' @param col_name Column to plot
#' @export
plot_density_per_column <- function(eval_list, col_name) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }

  # Build combined plotting data and drop methods missing the column
  pieces <- lapply(names(eval_list), function(method) {
    e <- eval_list[[method]]
    if (is.null(e$true[[col_name]]) || is.null(e$imputed[[col_name]])) return(NULL)
    n_true <- length(e$true[[col_name]])
    n_imp  <- length(e$imputed[[col_name]])
    if (n_true == 0 || n_imp == 0) return(NULL)

    data.frame(
      Value  = c(as.numeric(e$true[[col_name]]), as.numeric(e$imputed[[col_name]])),
      Type   = factor(c(rep("True", n_true), rep("Imputed", n_imp))),
      Method = factor(method),
      stringsAsFactors = FALSE
    )
  })

  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0) {
    stop("No data available for column '", col_name, "' across the provided methods.")
  }

  plot_data <- do.call(rbind, pieces)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Value, fill = .data$Type)) +
    ggplot2::geom_density(alpha = 0.4, colour = NA) +
    ggplot2::facet_wrap(~ .data$Method, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("True vs Imputed - Density for", col_name),
      x = col_name, y = "Density"
    ) +
    ggplot2::theme(legend.position = "top")
}

#' Plot consolidated density comparisons for all columns and methods
#' @param eval_list Output of get_eval_list(res)
#' @export
plot_density_all <- function(eval_list) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }

  # Stack everything: Method x Column x (True/Imputed)
  pieces <- lapply(names(eval_list), function(method) {
    e <- eval_list[[method]]
    cols <- intersect(names(e$true), names(e$imputed))
    if (length(cols) == 0) return(NULL)

    do.call(rbind, lapply(cols, function(col) {
      tv <- e$true[[col]]; iv <- e$imputed[[col]]
      tv <- as.numeric(tv); iv <- as.numeric(iv)
      tv <- tv[is.finite(tv)]; iv <- iv[is.finite(iv)]
      if (length(tv) == 0 || length(iv) == 0) return(NULL)

      data.frame(
        Value  = c(tv, iv),
        Type   = factor(c(rep("True", length(tv)), rep("Imputed", length(iv)))),
        Method = factor(method),
        Column = factor(col),
        stringsAsFactors = FALSE
      )
    }))
  })

  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0) {
    stop("No evaluable columns found to plot.")
  }

  all_data <- do.call(rbind, pieces)

  ggplot2::ggplot(all_data, ggplot2::aes(x = .data$Value, fill = .data$Type)) +
    ggplot2::geom_density(alpha = 0.35, colour = NA) +
    ggplot2::facet_grid(.data$Column ~ .data$Method, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "True vs Imputed - Density Comparisons (All Columns x Methods)",
      x = "Value", y = "Density"
    ) +
    ggplot2::theme(legend.position = "top")
}

