# imputetoolkit <a href="https://tanveer09.github.io/imputetoolkit/"><img src="man/figures/imputetoolkit_logo.png" align="right" height="120" width="120" style="margin-top:-5px;"/></a>


<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/tanveer09/imputetoolkit/blob/main/LICENSE.md)
[![Docs: pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://github.com/tanveer09/imputetoolkit)
[![GitHub repo](https://img.shields.io/badge/GitHub-imputetoolkit-black?logo=github)](https://github.com/tanveer09/imputetoolkit)
[![Made with Rcpp & R](https://img.shields.io/badge/Made%20with-Rcpp%20%26%20R-blue.svg)]()
<!-- badges: end -->



## Overview

**imputetoolkit** is an R package for **evaluating and benchmarking missing-data imputation techniques** using a unified **R + C++ backend**.

It automates the full workflow of **data preparation → missingness injection → imputation → metric-based evaluation → visualization**, providing both **numeric** and **categorical** assessments.

With parallel computation, C++ acceleration, and an intuitive S3 interface (`print`, `summary`, `plot`), the package allows users to **quantitatively justify imputation strategy choices**.

------------------------------------------------------------------------

## Key Features

| Category                      | Description                                                                                                                      |
|----------------------|-------------------------------------------------|
| **Imputation Methods**        | Mean/Mode (skewness-aware log transform), Median/Mode, MICE, and mixed-type KNN (parallelized)        |
| **Evaluation Metrics (Rcpp)** | Split into: <br>• *Numeric:* RMSE, MAE, R², Correlation, KS <br>• *Categorical:* Accuracy, Kappa, F1, MacroF1, Balanced Accuracy |
| **Parallelization**           | Multi-core KNN imputation using `FNN`, `foreach`, and `doParallel`                                                               |
| **Scaling**                   | Min-max scaling of numeric columns for fair metric comparison                                                                    |
| **Visualization**             | ggplot2-based metric and density plots (per method or all methods)                                                               |
| **Recommendation**            | Automatic “best-method” suggestions per metric type                                                                              |
| **Reproducibility**           | Deterministic pipeline with fixed seeds and consistent outputs                                                                   |

------------------------------------------------------------------------

## Package Structure

```         

imputetoolkit/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── evaluator.R             # main pipeline + S3 interface
│   ├── plot_metrics.R          # plotting utilities
│   ├── print_summary.R         # print() & summary() methods
│   ├── utils.R                 # helper utilities
│   └── RcppExports.R           # glue between R and C++
├── src/
│   ├── evaluate_imputation_split.cpp   # C++ backend for numeric + categorical metrics
│   ├── RcppExports.cpp
│   └── (future extensions)
├── inst/extdata/               # sample datasets
├── vignettes/                  # user guide (html_vignette)
└── tests/testthat/             # full unit test suite
```

------------------------------------------------------------------------

## Core Functions

### `evaluator()`

Main pipeline that executes: 1. Data loading or input validation\
2. Missingness injection\
3. Imputation via all four methods\
4. Metric evaluation (numeric & categorical)\
5. Return of structured `evaluator` objects

**Returns:**

``` r
list(mean_mode, median_mode, mice, knn)
```

Each object contains:

``` r
$method
$metrics_numeric
$metrics_categorical
```

### S3 Methods

-   `print.evaluator()` – concise summary by metric type
-   `summary.evaluator()` – per-column and aggregated results
-   `print_metrics()` – combined table across all methods
-   `plot_metrics()` – ggplot visualization for any metric
-   `suggest_best_method()` – identifies top-performing methods
-   `evaluate_results()` – runs print + plot + recommendation
-   `get_eval_list()` – extracts true/imputed pairs for density plots
-   `plot_density_per_column()` / `plot_density_all()` – visualize true vs imputed distributions

------------------------------------------------------------------------

## C++ Backend (`evaluate_imputation_split.cpp`)

Efficient Rcpp backend that calculates both **numeric** and **categorical** metrics in one pass.

### Numeric metrics

| Metric      | Description                  |
|:------------|:-----------------------------|
| RMSE        | Root Mean Square Error       |
| MAE         | Mean Absolute Error          |
| R²          | Coefficient of Determination |
| Correlation | Pearson correlation          |
| KS          | Kolmogorov–Smirnov statistic |

### Categorical metrics

| Metric           | Description               |
|:-----------------|:--------------------------|
| Accuracy         | Exact match ratio         |
| Kappa            | Chance-adjusted agreement |
| F1               | Micro-F1 score            |
| MacroF1          | Mean F1 over all classes  |
| BalancedAccuracy | Mean recall per class     |

All metrics are implemented in C++ using `Rcpp`, ensuring fast and stable performance.

------------------------------------------------------------------------

## Unit Tests (`tests/testthat/`)

Comprehensive suite covering:

| Test Area             | Focus                                                    |
|------------------------|-----------------------------------------------|
| **Pipeline**          | Verifies all 4 imputation methods run successfully       |
| **S3 Methods**        | `print`, `summary`, and invisibility checks              |
| **Wrapper Functions** | `extract_metrics`, `plot_metrics`, `suggest_best_method` |
| **Error Handling**    | Invalid inputs, unsupported file formats                 |
| **Reproducibility**   | Consistent metrics with fixed seed                       |
| **Parallel KNN**      | Multi-core execution validation                          |
| **C++ Evaluator**     | Direct unit tests for numeric/categorical metrics        |

Example:

```         
[ FAIL 0 | WARN 1 | SKIP 0 | PASS 76 ]
```

------------------------------------------------------------------------

## Installation

Requires **R ≥ 4.0**, Rtools (Windows), and compilation support for Rcpp.

### Install from GitHub

``` r
install.packages("remotes")
remotes::install_github("tanveer09/imputetoolkit@main", build_vignettes = TRUE, INSTALL_opts = c("--install-tests"))
```

### Load the package

``` r
library(imputetoolkit)
```

### Browse documentation

``` r
browseVignettes("imputetoolkit")
```

------------------------------------------------------------------------

## Example Usage

### 1. Load Data

``` r
file <- system.file("extdata", "synthetic_mixed_dataset.csv", package = "imputetoolkit")
raw_data <- read.csv(file, stringsAsFactors = TRUE)
```

### 2. Run the Evaluator

``` r
res <- evaluator(data = raw_data)
```

### 3. Inspect One Method

``` r
print(res$mean_mode)
summary(res$mean_mode)
```

### 4. Compare All Methods

``` r
print_metrics(res)
plot_metrics(res, "ALL")
```

### 5. Suggest the Best Method

``` r
suggest_best_method(res, "ALL")
```

Example output:

```         
Numeric Columns:
     Best imputation method as per "RMSE, MAE, R2, KS" metric: Mean/Mode
     Best imputation method as per "Correlation" metric: MICE

Categorical Columns:
     Best imputation method as per "Kappa" metric: KNN
     Best imputation method as per "Accuracy, F1, BalancedAccuracy" metric: Mean/Mode
     Best imputation method as per "MacroF1" metric: MICE
```

------------------------------------------------------------------------

## Example Visuals

### Metric Comparison Across Methods

``` r
plot_metrics(res, "ALL")
```

![Comparison of Imputation Methods across All Metrics](docs/plots/example_metrics.png)

### Single Metric (e.g. RMSE)

``` r
plot_metrics(res, "RMSE")
```

### Density Comparison

``` r
eval_list <- get_eval_list(res)
plot_density_per_column(eval_list, "age")
plot_density_all(eval_list)
```

![True vs Imputed Densities](docs/plots/example_density.png)

------------------------------------------------------------------------

## Evaluation Metrics Summary

| Metric               | Type        | Goal | Description                        |
|:-----------------|:-----------------|:----------------:|:------------------|
| **RMSE**             | Numeric     |  ↓   | Root Mean Squared Error            |
| **MAE**              | Numeric     |  ↓   | Mean Absolute Error                |
| **R²**               | Numeric     |  ↑   | Variance explained                 |
| **Correlation**      | Numeric     |  ↑   | Pearson correlation                |
| **KS**               | Numeric     |  ↑   | Distribution similarity            |
| **Accuracy**         | Categorical |  ↑   | Exact match rate                   |
| **Kappa**            | Categorical |  ↑   | Chance-corrected agreement         |
| **F1 / MacroF1**     | Categorical |  ↑   | Balance between precision & recall |
| **BalancedAccuracy** | Categorical |  ↑   | Mean recall per class              |

------------------------------------------------------------------------

## Additional Features

-   **Split-metric evaluation:** numeric and categorical metrics handled separately.
-   **C++ acceleration:** high-performance evaluation with safe NA handling.
-   **Parallel KNN:** multi-core computation for mixed-type datasets.
-   **Skewness-aware mean imputation:** applies log + geometric mean for \|skewness\| \> 1.
-   **Reproducible pipeline:** consistent metrics across runs (`set.seed()` used).
-   **Custom plotting:** side-by-side bar charts, density overlays, and panel comparisons.

------------------------------------------------------------------------

## Documentation

-   Function references: `?evaluator`, `?plot_metrics`, `?suggest_best_method`
-   Vignette tutorial: `vignette("imputetoolkit")`
-   Online docs: [pkgdown site](https://tanveer09.github.io/imputetoolkit/)

------------------------------------------------------------------------

## Development Notes

-   Modular design: easy extension for new imputation algorithms (e.g. missForest, EM).
-   Full `testthat` coverage with \>50 passing tests.
-   Fully vectorized Rcpp implementation.
-   Works seamlessly with `knitr` / `rmarkdown` for reproducible reports.

------------------------------------------------------------------------

## Feedback and Contributions

Feedback, feature requests, or bug reports are welcome via the [Issues page](https://github.com/tanveer09/imputetoolkit/issues).

------------------------------------------------------------------------

## License

Released under the [MIT License](LICENSE).

------------------------------------------------------------------------

## Citation

> Singh, Tanveer. (2025). *imputetoolkit: An R Package for Evaluating Missing Data Imputation Methods.* Victoria University of Wellington.

------------------------------------------------------------------------
