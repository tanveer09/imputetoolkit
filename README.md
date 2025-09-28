# imputetoolkit <img src="https://www.r-project.org/logo/Rlogo.png" width="40" align="right"/>

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![pkgdown site](https://img.shields.io/badge/docs-pkgdown-blue)](https://yourusername.github.io/imputetoolkit/) [![Made with R](https://img.shields.io/badge/Made%20with-Rcpp%20%26%20R-blue.svg)]()

| An R package for **evaluating missing data imputation methods** with a unified **R + C++ backend**. It provides a full pipeline for **data cleaning, missingness injection, imputation, and evaluation**, with support for **multiple imputation methods** and a suite of **evaluation metrics**. \|
| Implemented with an efficient **Rcpp backend** and an intuitive **S3 interface** (`print`, `summary`, `plot`), the package is designed for **benchmarking imputation strategies** and teaching/experimentation with missing data handling. \|

## Package Components

### Core Structure

-   **DESCRIPTION / NAMESPACE**: Standard package metadata, dependencies, and exported functions.
-   **R/**: R interface functions, wrapper utilities, and S3 methods.
-   **src/**: C++ backend (`evaluator.cpp`, `imputer.cpp`) compiled via Rcpp.
-   **tests/testthat/**: Unit tests for R and C++ integration.
-   **vignettes/**: Demonstration of package usage and design.
-   **inst/extdata/**: Sample dataset (`sample_dataset.csv`).

### R Files

-   **api.R**
    -   Entry point functions for users.
-   **imputer.R**
    -   Implements helper routines for imputation (mean, median, mode, MICE).
    -   **evaluator()**
        -   Main pipeline: `evaluator()`\
        -   Steps:
            -   Load dataset (CSV, TSV, TXT, Excel, RDS).\
            -   Clean categorical variables, inject controlled missingness.\
            -   Apply multiple imputation strategies: **Mean/Mode**, **Median/Mode**, **MICE**.\
            -   Collect true vs. imputed values.\
            -   Pass to the **C++ backend evaluator** for metric computation.\
        -   Exports S3 methods:
            -   `print.evaluator`: concise textual summary.\
            -   `summary.evaluator`: per-column + global metrics table.\
            -   `plot_metrics`: ggplot comparison of methods.\
            -   `print_metrics`: table-formatted comparison via knitr.\
            -   `suggest_best_method`: choose method by metric (min/max).\
            -   `evaluate_results`: wrapper to print, plot, and suggest best method.
-   **RcppExports.R**
    -   Auto-generated glue between R and C++ (via Rcpp).

### C++ Files (src/)

-   **evaluator.cpp**
    -   Implements `evaluate_imputation()` (exported to R).\
    -   Defines an **Evaluator class** that computes:
        -   RMSE, MAE, R², correlation, KS statistic, Accuracy.\
        -   Per-column metrics and global averages.\
    -   Uses **OOP encapsulation** (private methods for column-level evaluation, public interface returning results).
-   **imputer.cpp** *(extension point)*
    -   Placeholder for additional algorithms (future work).
-   **RcppExports.cpp**
    -   Auto-generated bindings between R and C++.


### Testing (tests/testthat/)

-   **test-evaluator.R**
    -   Tests the full pipeline with the synthetic dataset.\
    -   Verifies all methods (`mean_mode`, `median_mode`, `mice`).\
    -   Checks:
        -   Metrics existence (`RMSE`, `MAE`, `R2`, `Correlation`, `KS`, `Accuracy`).\
        -   `print` and `summary` invisibility.\
        -   Wrapper functions (`extract_metrics`, `print_metrics`, `plot_metrics`, `suggest_best_method`).\
        -   Error handling (missing filename/data, unsupported file types, invalid metric names).\
        -   Reproducibility with fixed seed.
-   **test-rcpp-evaluator.R**
    -   Tests `evaluate_imputation()` directly with toy numeric data.\
    -   Validates metrics computation, correct method labeling, and output structure.

Together, these tests ensure **robustness, reproducibility, and informative error handling**.

------------------------------------------------------------------------

## Installation

Make sure you have R (≥ 4.0) and the **devtools/remotes** package:

``` r
install.packages("devtools")
devtools::install_github("tanveer09/imputetoolkit@draft")
```

Include unit tests during installation and execute them:

``` r
install.packages("remotes")
remotes::install_github("tanveer09/imputetoolkit@draft", INSTALL_opts = c("--install-tests"), force = TRUE)

library(testthat)
test_package("imputetoolkit")
```


Then load:

``` r
library(imputetoolkit)
```

------------------------------------------------------------------------

## Sample Dataset

The package ships with a **synthetic mixed-type dataset** located at:

```         
inst/extdata/sample_dataset.csv
```

This dataset contains **numeric and categorical columns** with controlled missingness, useful for benchmarking imputation strategies.

### Example

``` r
filename <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
raw_data <- read.csv(filename, stringsAsFactors = TRUE)
res <- evaluator(data = raw_data)
summary(res$mean_mode)
#metrics_df <- extract_metrics(res)
print_metrics(res)
plot_metrics(res, "RMSE")
suggest_best_method(res, "Accuracy", higher_better = TRUE)
```

------------------------------------------------------------------------

## Usage Examples

### Run Evaluator

``` r
res <- evaluator(filename = file)
```

### Print & Summarize

``` r
print(res$mean_mode)
summary(res$mean_mode)
```

### Compare Methods

``` r
print_metrics(res)
plot_metrics(res, "RMSE")
```

### Suggest Best Method

``` r
suggest_best_method(res, metric = "RMSE", higher_better = FALSE)
```

------------------------------------------------------------------------

## Documentation

-   Function manuals: `?evaluator`, `?plot_metrics`, etc.
-   Vignettes: under `vignettes/` (in progress).
-   HTML docs: pkgdown site at <https://tanveer09.github.io/imputetoolkit/>

------------------------------------------------------------------------

## Development Notes

-   C++ backend (`Evaluator` class) demonstrates **OOP encapsulation**.
-   S3 interface on the R side for consistent `print`, `summary`, and plotting.
-   Roadmap includes adding **kNN** and **missForest** imputations.
-   Current release: **core pipeline stable + tested**.

------------------------------------------------------------------------

## License

This package is released under the [MIT License](LICENSE).

------------------------------------------------------------------------
