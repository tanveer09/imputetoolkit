#' evaluator
#'
#' Full pipeline: data cleaning, missingness injection, imputation, and evaluation.
#'
#' @param data is either a data.frame that has already been loaded in R. or is the filename to a data file (CSV, TSV, TXT, XLSX, XLS, RDS)
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
#' @importFrom stats sd

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
evaluator <- function(data = NULL) {
  set.seed(123) # reproducibility

  ##########################
  # ---- 0. Load data ---- #
  ##########################

  raw_data <- data

  if (is.character(data)) {
    file_ext <- tools::file_ext(data)

    if (file_ext == "csv") {
      raw_data <- utils::read.csv(data, stringsAsFactors = TRUE)
    } else if (file_ext %in% c("tsv", "txt")) {
      raw_data <- utils::read.delim(data, stringsAsFactors = TRUE)
    } else if (file_ext %in% c("xlsx", "xls")) {
      raw_data <- readxl::read_excel(data)
    } else if (file_ext %in% c("rds")) {
      raw_data <- readRDS(data)
    } else {
      stop("Unsupported file type: ", file_ext)
    }
  }

  if (is.null(data)) {
    stop("No input provided. Please pass either:
  - a data.frame (e.g., evaluator(mydata)), or
  - a file path (e.g., evaluator('data.csv')).")
  }

  if (!is.data.frame(raw_data)) {
    stop("Input must resolve to a data.frame after loading.")
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

  # ---- Stop if no missing values found ----
  if (length(raw_data_col_missingness) == 0) {
    stop("No missing values detected in the dataset.
       Please provide a dataset with at least one missing value for imputation.")
  }

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

  data_mice <- NULL
  if (ncol(raw_data_modified) >= 2) {
    tryCatch({
      mice_imp <- mice::mice(raw_data_modified, m = 1, method = "pmm", maxit = 5, seed = 123, printFlag = FALSE)
      data_mice <- mice::complete(mice_imp)
    }, error = function(e) {
      warning("MICE imputation skipped: ", conditionMessage(e))
      # fallback: use median/mode as a safe default
      data_mice <- data_median
    })
  } else {
    warning("MICE skipped because dataset has fewer than 2 columns.")
    data_mice <- data_median
  }

  ##############################################
  # ---- KNN Imputation (robust version) ---- #
  ##############################################
  message("Running mixed-type parallel KNN imputation using FNN...")

  if (!requireNamespace("FNN", quietly = TRUE))
    stop("Package 'FNN' required for fast KNN imputation. Please install it.")
  if (!requireNamespace("doParallel", quietly = TRUE))
    stop("Package 'doParallel' required for parallel execution. Please install it.")
  if (!requireNamespace("foreach", quietly = TRUE))
    stop("Package 'foreach' required.")

  data_knn <- raw_data_modified

  num_cols <- names(which(sapply(data_knn, is.numeric)))
  cat_cols <- names(which(!sapply(data_knn, is.numeric)))

  # --- Standardize numeric data safely ---
  num_data <- data_knn[, num_cols, drop = FALSE]
  num_data_scaled <- scale(num_data)
  pred_data <- as.data.frame(num_data_scaled)

  # --- Setup parallel backend (safe) ---
  total_cores <- parallel::detectCores(logical = TRUE)
  is_check_env <- isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE"))) ||
    ("CHECK" %in% toupper(names(Sys.getenv()))) ||
    ("TESTTHAT" %in% toupper(names(Sys.getenv())))
  cores <- if (is_check_env) 1 else max(1, total_cores - 1)

  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    message("Using ", cores, " cores for KNN imputation...")
  } else {
    foreach::registerDoSEQ()
    message("Using sequential mode for KNN imputation (testing environment)...")
  }

  safe_knn <- function(train_x, test_x, k = 3) {
    if (nrow(train_x) == 0 || nrow(test_x) == 0)
      return(list(nn.index = matrix(integer(0), nrow = 0, ncol = k)))
    tryCatch(
      FNN::get.knnx(train_x, test_x, k = k),
      error = function(e) list(nn.index = matrix(integer(0), nrow = 0, ncol = k))
    )
  }

  # ---- Impute numeric columns ----
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
    test_x  <- pred_data[test_idx,  pred_cols, drop = FALSE]

    # Skip if predictors have zero variance
    if (any(apply(train_x, 2, stats::sd, na.rm = TRUE) == 0)) return(col_data)

    knn_res <- safe_knn(train_x, test_x, k = 3)
    if (nrow(knn_res$nn.index) == 0) return(col_data)

    imp_vals <- vapply(seq_len(nrow(knn_res$nn.index)), function(i) {
      neighbor_idx <- knn_res$nn.index[i, ]
      mean(col_data[train_idx[neighbor_idx]], na.rm = TRUE)
    }, numeric(1))

    # replace only missing entries
    col_data[test_idx] <- imp_vals
    return(col_data)
  }

  # ---- Reconstruct numeric data from scaled form ----
  if (length(num_cols) > 0 && length(imputed_numeric) > 0) {
    # Convert back to a named data.frame to avoid lost column names
    num_data_scaled_df <- as.data.frame(num_data_scaled)
    for (i in seq_along(num_cols)) {
      colname <- num_cols[i]
      vec <- imputed_numeric[[i]]
      if (!is.null(vec)) {
        # replace column values directly even if some NAs remain
        num_data_scaled_df[[colname]] <- vec
      }
    }

    # Retrieve center/scale attributes
    attr_center <- attr(num_data_scaled, "scaled:center")
    attr_scale  <- attr(num_data_scaled, "scaled:scale")

    # Reverse scaling correctly (column-by-column)
    num_data_imputed <- as.data.frame(lapply(seq_along(num_cols), function(i) {
      colname <- num_cols[i]
      v <- num_data_scaled_df[[colname]]
      sc <- attr_scale[colname]
      ct <- attr_center[colname]
      if (is.null(sc) || is.na(sc) || sc == 0) return(v)
      return(v * sc + ct)
    }))
    names(num_data_imputed) <- num_cols

    # Write back to main data
    for (colname in num_cols) {
      data_knn[[colname]] <- num_data_imputed[[colname]]
    }
  }

  cat("\n[DEBUG] KNN numeric NAs before:", sum(is.na(raw_data_modified[num_cols])),
      "after:", sum(is.na(data_knn[num_cols])), "\n")


  # ---- Impute categorical columns ----
  imputed_categorical <- foreach::foreach(col = cat_cols, .packages = "FNN") %dopar% {
    col_data <- data_knn[[col]]
    na_idx <- which(is.na(col_data))
    if (length(na_idx) == 0) return(as.character(col_data))

    if (length(num_cols) == 0) return(as.character(col_data))

    train_idx <- which(!is.na(col_data) & complete.cases(pred_data))
    test_idx  <- which(is.na(col_data) & complete.cases(pred_data))
    if (length(train_idx) == 0 || length(test_idx) == 0)
      return(as.character(col_data))

    train_x <- pred_data[train_idx, , drop = FALSE]
    test_x  <- pred_data[test_idx,  , drop = FALSE]

    if (any(apply(train_x, 2, sd, na.rm = TRUE) == 0)) return(as.character(col_data))

    knn_res <- safe_knn(train_x, test_x, k = 3)
    if (nrow(knn_res$nn.index) == 0) return(as.character(col_data))

    get_mode <- function(x) {
      ux <- unique(x[!is.na(x)])
      if (length(ux) == 0) return(NA_character_)
      ux[which.max(tabulate(match(x, ux)))]
    }

    imp_vals <- vapply(seq_len(nrow(knn_res$nn.index)), function(i) {
      neighbor_idx <- knn_res$nn.index[i, ]
      get_mode(as.character(col_data[train_idx[neighbor_idx]]))
    }, character(1))

    col_data[test_idx] <- imp_vals
    as.character(col_data)
  }

  # ---- Merge categorical results back safely ----
  if (length(cat_cols) > 0 && length(imputed_categorical) > 0) {
    for (i in seq_along(cat_cols)) {
      colname <- cat_cols[i]
      imp_vals <- imputed_categorical[[i]]
      if (!is.null(imp_vals) && length(imp_vals) == nrow(data_knn)) {
        data_knn[[colname]] <- factor(imp_vals, levels = levels(raw_data[[colname]]))
      }
    }
  }

  # ---- Stop cluster safely ----
  if (exists("cl") && inherits(cl, "cluster")) {
    parallel::stopCluster(cl)
    message("Parallel mixed-type KNN imputation complete (parallel).")
  } else {
    message("KNN imputation completed in sequential mode (no cluster to stop).")
  }

  cat("\n[DEBUG] KNN numeric NA count before:", sum(is.na(raw_data_modified[num_cols])),
      "after:", sum(is.na(data_knn[num_cols])), "\n")



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


  ###############################################################
  # ---- 5. Call Rcpp evaluator (split numeric vs categorical) ----
  ###############################################################
  numeric_cols <- names(Filter(is.numeric, raw_data))
  categorical_cols <- names(Filter(Negate(is.numeric), raw_data))

  # ---- 5. Evaluate imputation (split numeric vs categorical) ----
  res_mean <- evaluate_imputation_split(eval_mean$true, eval_mean$imputed,
                                        numeric_cols, categorical_cols, "Mean/Mode")
  res_median <- evaluate_imputation_split(eval_median$true, eval_median$imputed,
                                          numeric_cols, categorical_cols, "Median/Mode")
  res_mice <- evaluate_imputation_split(eval_mice$true, eval_mice$imputed,
                                        numeric_cols, categorical_cols, "MICE")
  res_knn <- evaluate_imputation_split(eval_knn$true, eval_knn$imputed,
                                       numeric_cols, categorical_cols, "KNN")

  # ---- Attach true vs imputed evaluation data for density plots ----
  res_mean$eval   <- eval_mean
  res_median$eval <- eval_median
  res_mice$eval   <- eval_mice
  res_knn$eval    <- eval_knn

  # assign class
  class(res_mean) <- class(res_median) <- class(res_mice) <- class(res_knn) <- "evaluator"

  # ---- Return complete result ----
  list(
    mean_mode   = res_mean,
    median_mode = res_median,
    mice        = res_mice,
    knn         = res_knn
  )
}

#' Print evaluator results separately for numeric and categorical
#'
#' @title Print Evaluator Results
#' @description Prints numeric and categorical evaluation metrics for an evaluator object.
#'
#' @param x An object of class \code{evaluator}, typically returned by \code{evaluator()}.
#' @param ... Additional arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns the evaluator object.
#' @method print evaluator
#' @export
print.evaluator <- function(x, ...) {
  cat("Evaluation for method:", x$method, "\n")

  cat("\n--- Numeric Metrics ---\n")
  if (length(x$metrics_numeric)) {
    num_tbl <- do.call(rbind, lapply(x$metrics_numeric, as.data.frame))
    print(round(colMeans(num_tbl, na.rm = TRUE), 4))
  } else {
    cat("(none)\n")
  }

  cat("\n--- Categorical Metrics ---\n")
  if (length(x$metrics_categorical)) {
    cat_tbl <- do.call(rbind, lapply(x$metrics_categorical, as.data.frame))
    print(round(colMeans(cat_tbl, na.rm = TRUE), 4))
  } else {
    cat("(none)\n")
  }

  invisible(x)
}


#' @title Summarize Evaluator Results
#' @description Returns numeric and categorical metric tables for an evaluator object.
#'
#' @param object An object of class \code{evaluator}, typically returned by \code{evaluator()}.
#' @param ... Additional arguments passed to or from other methods (ignored).
#'
#' @return A list containing:
#'   \item{numeric_metrics}{A data frame of numeric metric values per column.}
#'   \item{categorical_metrics}{A data frame of categorical metric values per column.}
#' @method summary evaluator
#' @export
summary.evaluator <- function(object, ...) {
  if (length(object$metrics_numeric)) {
    num_df <- do.call(rbind, lapply(object$metrics_numeric, as.data.frame))
    num_df <- cbind(Column = names(object$metrics_numeric), num_df)
    rownames(num_df) <- NULL
  } else {
    num_df <- data.frame()
  }

  if (length(object$metrics_categorical)) {
    cat_df <- do.call(rbind, lapply(object$metrics_categorical, as.data.frame))
    cat_df <- cbind(Column = names(object$metrics_categorical), cat_df)
    rownames(cat_df) <- NULL
  } else {
    cat_df <- data.frame()
  }

  list(
    numeric_metrics = num_df,
    categorical_metrics = cat_df
  )
}


#' extract_metrics
#' Extract evaluation metrics (handles split numeric/categorical)
#' @keywords internal
#' @noRd
extract_metrics <- function(res) {
  if (!is.list(res)) stop("res must be a list returned by evaluator().")

  metrics_list <- lapply(names(res), function(method_name) {
    method_res <- res[[method_name]]

    # --- Aggregate numeric metrics ---
    if (!is.null(method_res$metrics_numeric) && length(method_res$metrics_numeric) > 0) {
      num_tbl <- do.call(rbind, lapply(method_res$metrics_numeric, as.data.frame))
      num_summary <- colMeans(num_tbl, na.rm = TRUE)
    } else {
      num_summary <- c()
    }

    # --- Aggregate categorical metrics ---
    if (!is.null(method_res$metrics_categorical) && length(method_res$metrics_categorical) > 0) {
      cat_tbl <- do.call(rbind, lapply(method_res$metrics_categorical, as.data.frame))
      cat_summary <- colMeans(cat_tbl, na.rm = TRUE)
    } else {
      cat_summary <- c()
    }

    # Combine numeric + categorical summaries
    all_metrics <- c(num_summary, cat_summary)

    data.frame(
      Method = method_res$method,
      t(as.data.frame(all_metrics)),
      stringsAsFactors = FALSE
    )
  })

  metrics_df <- do.call(rbind, metrics_list)
  rownames(metrics_df) <- NULL
  metrics_df
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
#' Plot evaluation metrics (numeric + categorical)
#'
#' @param x A result list returned by evaluator()
#' @param metric Character, which metric to plot (e.g., "RMSE", "MAE", "R2", "Accuracy", "F1", "ALL")
#' @return A ggplot2 object
#' @export
plot_metrics <- function(x, metric = "RMSE") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required. Please install it.")
  }

  # ---- Extract all available metrics ----
  metrics_df <- extract_metrics(x)

  # ---- Detect numeric + categorical metric columns ----
  numeric_metrics <- c("RMSE", "MAE", "R2", "Correlation", "KS")
  categorical_metrics <- c("Accuracy", "Kappa", "F1",
                           "MacroF1", "BalancedAccuracy")

  all_metrics <- intersect(
    c(numeric_metrics, categorical_metrics),
    colnames(metrics_df)
  )

  if (length(all_metrics) == 0) {
    stop("No valid metrics found in the results.")
  }

  # ---- When user wants all metrics ----
  if (metric == "ALL") {
    long_df <- tidyr::pivot_longer(
      metrics_df,
      cols = dplyr::all_of(all_metrics),
      names_to = "Metric",
      values_to = "Value"
    )

    p <- ggplot2::ggplot(
      long_df,
      ggplot2::aes(x = .data$Method, y = .data$Value, fill = .data$Method)
    ) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::facet_wrap(~ Metric, scales = "free_y") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Comparison of Imputation Methods across All Metrics",
        x = "Method",
        y = "Value"
      ) +
      ggplot2::theme(legend.position = "none")

    return(p)
  }

  # ---- When user specifies a single metric ----
  if (!(metric %in% all_metrics)) {
    stop(
      "Metric '", metric, "' not found. Available metrics are: ",
      paste(all_metrics, collapse = ", "), " or use metric = 'ALL'."
    )
  }

  p <- ggplot2::ggplot(
    metrics_df,
    ggplot2::aes(x = .data$Method, y = .data[[metric]], fill = .data$Method)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Comparison of Methods on", metric),
      x = "Method",
      y = metric
    ) +
    ggplot2::theme(legend.position = "none")

  return(p)
}


#' Suggest the best imputation method separately for numeric and categorical columns
#'
#' @param x A result list returned by evaluator()
#' @param metric Character, metric to optimize
#'   ("RMSE", "MAE", "R2", "Correlation", "KS", "Accuracy", "Kappa",
#'    "F1", "MacroF1", "BalancedAccuracy", or "ALL").
#'   If "ALL" or missing, compares across all metrics.
#' @return A list with two elements: `numeric` and `categorical`,
#'   each containing the best methods based on the chosen metric(s).
#' @export
suggest_best_method <- function(x, metric = "ALL") {
  metrics_df <- extract_metrics(x)

  # ---- Define metric groups ----
  numeric_metrics <- c("RMSE", "MAE", "R2", "Correlation", "KS")
  categorical_metrics <- c("Accuracy", "Kappa", "F1", "MacroF1", "BalancedAccuracy")

  higher_better <- c("R2", "Correlation", "KS", "Accuracy", "Kappa", "F1", "MacroF1", "BalancedAccuracy")
  lower_better  <- c("RMSE", "MAE")

  # ---- Helper to select best method ----
  pick_best <- function(df, metric) {
    if (!metric %in% colnames(df)) return(NA_character_)
    vals <- df[[metric]]
    if (all(is.na(vals))) return(NA_character_)
    if (metric %in% higher_better) {
      df$Method[which.max(vals)]
    } else {
      df$Method[which.min(vals)]
    }
  }

  # ---- Split numeric vs categorical ----
  num_df <- metrics_df[, c("Method", intersect(numeric_metrics, colnames(metrics_df))), drop = FALSE]
  cat_df <- metrics_df[, c("Method", intersect(categorical_metrics, colnames(metrics_df))), drop = FALSE]

  # ---- Case 1: ALL metrics ----
  if (metric == "ALL" || is.null(metric)) {
    best_numeric <- lapply(intersect(numeric_metrics, colnames(num_df)), function(m) pick_best(num_df, m))
    best_categorical <- lapply(intersect(categorical_metrics, colnames(cat_df)), function(m) pick_best(cat_df, m))

    names(best_numeric) <- intersect(numeric_metrics, colnames(num_df))
    names(best_categorical) <- intersect(categorical_metrics, colnames(cat_df))

    # Group by method
    grouped_numeric <- split(names(best_numeric), unlist(best_numeric))
    grouped_categorical <- split(names(best_categorical), unlist(best_categorical))

    # ---- Nicely formatted output ----
    message("\nSuggested Best Imputation Methods:\n")

    if (length(grouped_numeric)) {
      message("Numeric Columns:")
      num_msg <- paste(vapply(names(grouped_numeric), function(m) {
        paste0("     Best imputation method as per \"",
               paste(grouped_numeric[[m]], collapse = ", "),
               "\" metric: ", m)
      }, character(1L)), collapse = "\n")
      message(num_msg)
    }

    if (length(grouped_categorical)) {
      message("\nCategorical Columns:")
      cat_msg <- paste(vapply(names(grouped_categorical), function(m) {
        paste0("     Best imputation method as per \"",
               paste(grouped_categorical[[m]], collapse = ", "),
               "\" metric: ", m)
      }, character(1L)), collapse = "\n")
      message(cat_msg)
    }

    return(invisible(list(
      numeric = grouped_numeric,
      categorical = grouped_categorical
    )))
  }

  # ---- Case 2: Single metric ----
  metric <- match.arg(metric, choices = c(numeric_metrics, categorical_metrics), several.ok = FALSE)

  if (metric %in% numeric_metrics) {
    best_method <- pick_best(num_df, metric)
    message("\nNumeric Columns:")
    message("     Best imputation method as per \"", metric, "\" metric: ", best_method)
    return(invisible(list(numeric = best_method, categorical = NULL)))
  } else if (metric %in% categorical_metrics) {
    best_method <- pick_best(cat_df, metric)
    message("\nCategorical Columns:")
    message("     Best imputation method as per \"", metric, "\" metric: ", best_method)
    return(invisible(list(numeric = NULL, categorical = best_method)))
  } else {
    stop("Unknown metric: ", metric)
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

