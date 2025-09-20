# raw_data_missing_row_indexes <- which(apply(is.na(raw_data), 1, any))
# length(raw_data_missing_row_indexes)
# 
# clean_data <- na.omit(raw_data)
# 
# no_missing_data <- raw_data[-raw_data_missing_row_indexes, ]
# 
# any_missingness <- colSums(is.na(no_missing_data)) / nrow(no_missing_data) * 100
# any_missingness <- any_missingness[col_missingness > 0]
# 

# 
# unique_missing_indices <- unique(unlist(missing_indices_per_col))
# length(unique_missing_indices)
# 
# total_missing <- sum(is.na(raw_data))
# total_missing
# 
# 
# 
# filename <- "synthetic_dataset.csv"
# 
# raw_data <- read.csv(filename)
# 
# raw_data_col_missingness <- colSums(is.na(raw_data)) / nrow(raw_data) * 100
# raw_data_col_missingness <- raw_data_col_missingness[raw_data_col_missingness > 0]
# 
# 
# raw_data_missing_indices_per_col <- lapply(raw_data, function(x) which(is.na(x)))
# raw_data_missing_indices_per_col



set.seed(123)  # for reproducibility

filename <- "synthetic_dataset.csv"
raw_data <- read.csv(filename)
dim(raw_data)
# 1. Get missingness % per column
raw_data_col_missingness <- colSums(is.na(raw_data)) / nrow(raw_data) * 100
raw_data_col_missingness <- raw_data_col_missingness[raw_data_col_missingness > 0]

# 2. For each column with missingness, get the sampled indices
raw_data_selected_indices <- lapply(names(raw_data_col_missingness), function(col) {
  
  # Get indices that are NOT missing
  non_missing_idx <- which(!is.na(raw_data[[col]]))
  
  # % missing for this column
  perc <- raw_data_col_missingness[col]
  
  # Number of indices to sample
  n_to_sample <- ceiling(length(non_missing_idx) * perc / 100)
  
  # Randomly sample
  sample(non_missing_idx, n_to_sample)
})

# Name the list elements by column name
names(raw_data_selected_indices) <- names(raw_data_col_missingness)

# 3. Make a copy and inject missing values
raw_data_modified <- raw_data

for (col in names(raw_data_selected_indices)) {
  idx <- raw_data_selected_indices[[col]]
  raw_data_modified[idx, col] <- NA
}

# raw_data_modified now has additional missing values
head(raw_data_modified)

modified_raw_data_col_missingness <- colSums(is.na(raw_data_modified)) / nrow(raw_data_modified) * 100
modified_raw_data_col_missingness <- modified_raw_data_col_missingness[modified_raw_data_col_missingness > 0]


# Function to compute mode
get_mode <- function(x) {
  ux <- unique(x[!is.na(x)])     # ignore NA
  ux[which.max(tabulate(match(x, ux)))]
}

# Copy dataset for imputation
#raw_data_imputed <- raw_data_modified

for (col in names(raw_data_modified)) {
  if (is.numeric(raw_data_modified[[col]])) {
    # Mean imputation for numeric
    mean_val <- mean(raw_data_modified[[col]], na.rm = TRUE)
    raw_data_modified[[col]][is.na(raw_data_modified[[col]])] <- mean_val
  } else {
    # Mode imputation for non-numeric
    mode_val <- get_mode(raw_data_modified[[col]])
    raw_data_modified[[col]][is.na(raw_data_modified[[col]])] <- mode_val
  }
}

# raw_data_imputed now has no missing values
#head(raw_data_modified)

#modified_raw_data_col_missingness <- colSums(is.na(raw_data_modified)) / nrow(raw_data_modified) * 100


# ---- Helper Functions ----
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

mae <- function(actual, predicted) {
  mean(abs(actual - predicted), na.rm = TRUE)
}

# ---- Metrics Calculation ----
metrics_per_col <- lapply(names(raw_data_selected_indices), function(col) {
  idx <- raw_data_selected_indices[[col]]
  actual <- raw_data[idx, col]
  imputed <- raw_data_modified[idx, col]
  
  if (is.numeric(actual)) {
    # For numeric columns
    data.frame(
      Column = col,
      RMSE = rmse(actual, imputed),
      MAE = mae(actual, imputed),
      Correlation = suppressWarnings(cor(actual, imputed, use = "complete.obs")),
      stringsAsFactors = FALSE
    )
  } else {
    # For categorical columns: only accuracy-like metrics
    accuracy <- mean(actual == imputed, na.rm = TRUE)
    data.frame(
      Column = col,
      RMSE = NA,
      MAE = NA,
      Correlation = accuracy,   # interpret correlation as %match
      stringsAsFactors = FALSE
    )
  }
})

# Combine into one data frame
metrics_df <- do.call(rbind, metrics_per_col)

# metrics_df contains per-column results
metrics_df
