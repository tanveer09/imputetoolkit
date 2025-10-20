# imputetoolkit <img src="https://www.r-project.org/logo/Rlogo.png" width="40" align="right"/>

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![pkgdown site](https://img.shields.io/badge/docs-pkgdown-blue)](https://tanveer09.github.io/imputetoolkit/)
[![Made with R](https://img.shields.io/badge/Made%20with-Rcpp%20%26%20R-blue.svg)]()

---

## Overview

**imputetoolkit** is an R package for **evaluating and benchmarking missing-data imputation methods** using a unified **R + C++ backend**.

It provides a complete, reproducible pipeline for:

- Data loading and cleaning  
- Controlled missingness injection  
- Multiple imputation techniques  
- Automated evaluation via Rcpp metrics  
- Visualization and comparison of results  

With parallel computation support and a clean S3 interface (`print`, `summary`, `plot`), **imputetoolkit** helps data scientists *quantitatively justify* their imputation strategy choices.

---

## Key Features

| Category | Description |
|-----------|-------------|
| **Imputation Methods** | Mean/Mode (with skewness-aware log transform), Median/Mode, MICE (Predictive Mean Matching), KNN (mixed-type, parallel) |
| **Evaluation Metrics (Rcpp)** | RMSE, MAE, R², Correlation, KS Statistic, Accuracy |
| **Scaling & Validation** | Min–max scaling of numeric columns before metric computation |
| **Visualization** | ggplot-based bar and density plots for all methods |
| **Recommendation** | Automatic best-method suggestion by metric |
| **Performance** | Parallelized KNN using `FNN` + `foreach` + `doParallel` |
| **Reproducibility** | Fixed seeds and consistent deterministic pipeline |

---

## Package Structure

```

imputetoolkit/
├── DESCRIPTION / NAMESPACE
├── R/
│   ├── evaluator.R          # main pipeline & helper functions
│   ├── plot_metrics.R       # plotting utilities
│   ├── print_summary.R      # S3 print/summary methods
│   ├── utils.R              # internal helpers
│   └── RcppExports.R        # auto-generated glue
├── src/
│   ├── evaluator.cpp        # Rcpp backend for metric computation
│   ├── RcppExports.cpp
│   └── (future extensions)
├── inst/extdata/            # sample datasets
├── vignettes/               # user guide (html_vignette)
└── tests/testthat/          # full unit test suite

```

---

## Core Functions

### `evaluator()`
Main entry point: runs the entire pipeline  
→ load data → inject missingness → impute → evaluate.

**Supported input types:** `.csv`, `.tsv`, `.xlsx`, `.rds`

**Returns:**  
A named list of `evaluator` S3 objects:  
`mean_mode`, `median_mode`, `mice`, `knn`

### S3 Methods
- `print.evaluator()` – concise textual summary  
- `summary.evaluator()` – per-column + global metrics  
- `print_metrics()` – tabular comparison via **knitr::kable**  
- `plot_metrics()` – ggplot visual comparison (single or all metrics)  
- `suggest_best_method()` – recommends optimal methods by metric  
- `evaluate_results()` – print + plot + recommend in one call  
- `get_eval_list()` – extract true vs imputed pairs for plotting  
- `plot_density_per_column()` / `plot_density_all()` – density overlays for distribution comparison  

---

## C++ Backend (`src/evaluator.cpp`)

Implements `evaluate_imputation()` via an OOP **Evaluator** class that:

- Computes per-column and global metrics:  
  RMSE, MAE, R², Correlation, KS, Accuracy  
- Performs safe mean aggregation (ignoring NAs)  
- Returns structured `List` objects to R

All numeric operations are optimized with **Rcpp**.

---

## Unit Tests (`tests/testthat/`)

Comprehensive tests verify:

| Test Area | Coverage |
|------------|-----------|
| **Core pipeline** | Mean/Mode, Median/Mode, MICE, and KNN |
| **S3 methods** | `print`, `summary`, invisibility checks |
| **Wrappers** | `extract_metrics`, `print_metrics`, `plot_metrics`, `suggest_best_method` |
| **Error handling** | Invalid inputs, unsupported file types, missing arguments |
| **Reproducibility** | Identical metrics with fixed seed |
| **Parallel KNN** | Sequential fallback in test environments |
| **Rcpp evaluator** | Independent numeric tests for correctness |

Example output:

```

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 46 ]

````

---

## Installation

Requires **R ≥ 4.0** and **Rtools** (Windows).

### Install via `remotes`
```r
install.packages("remotes")
remotes::install_github("tanveer09/imputetoolkit@draft", build_vignettes = TRUE)
````

### Load

```r
library(imputetoolkit)
```

### View Vignette

```r
browseVignettes("imputetoolkit")
```

---

## Example Usage

### 1. Load Data

```r
file <- system.file("extdata", "sample_dataset.csv", package = "imputetoolkit")
raw_data <- read.csv(file, stringsAsFactors = TRUE)
```

### 2. Run Evaluator

```r
res <- evaluator(data = raw_data)
```

### 3. Inspect Results

```r
print(res$mean_mode)
summary(res$mean_mode)
```

### 4. Compare All Methods

```r
print_metrics(res)
plot_metrics(res, "ALL")
```

### 5. Recommend Best Method

```r
suggest_best_method(res, "ALL")
```

### 6. Visualize True vs Imputed Densities

```r
eval_list <- get_eval_list(res)
plot_density_per_column(eval_list, "age")
plot_density_all(eval_list)
```

---

## Sample Output

```
Evaluation for method: KNN 
Global Metrics:
  RMSE       : 0.1841 
  MAE        : 0.1412 
  R^2        : 0.9213 
  Correlation: 0.9584 
  KS         : 0.8846 
  Accuracy   : 0.8125 

Per-column metrics available in x$metrics
```

---

## Evaluation Metrics

| Metric          | Description                                              | Goal |
| :-------------- | :------------------------------------------------------- | :--: |
| **RMSE**        | Root Mean Squared Error                                  |   ↓  |
| **MAE**         | Mean Absolute Error                                      |   ↓  |
| **R²**          | Proportion of variance explained                         |   ↑  |
| **Correlation** | Pearson correlation (true vs imputed)                    |   ↑  |
| **KS**          | Kolmogorov–Smirnov statistic (distributional similarity) |   ↑  |
| **Accuracy**    | Exact match rate (categorical)                           |   ↑  |

---

## Additional Features

* **Log-transform mean correction** for highly skewed numeric columns (`|skewness| > 1`)
* **Min–max normalization** before metric comparison
* **Parallel KNN** using all but one core (`FNN`, `foreach`, `doParallel`)
* **Sequential fallback** for CRAN/test environments
* **Density comparison plots** for visual validation

---

## Example Visuals

**Barplot of RMSE Across Methods**

```r
plot_metrics(res, "RMSE")
```

**Faceted Metrics Panel**

```r
plot_metrics(res, "ALL")
```

**Density Comparison Example**

```r
plot_density_per_column(eval_list, "income")
```

---

## Documentation

* Function references:
  `?evaluator`, `?plot_metrics`, `?suggest_best_method`
* Vignettes:
  `vignette("imputetoolkit")`
* Tutorials and examples:
  [pkgdown site](https://tanveer09.github.io/imputetoolkit/)

---

## Development Notes

* Fully **vectorized Rcpp** backend
* S3-based front-end for clean extensibility
* Unit-tested for reproducibility and accuracy
* Modular design: new imputation algorithms can be added easily
* Future extensions: missForest, EM, and deep autoencoder imputations

---

## License

Released under the [MIT License](LICENSE).

---

## Citation

> Singh, Tanveer. (2025). *imputetoolkit: An R Package for Evaluating Missing-Data Imputation Methods*.
> Victoria University of Wellington.

---
