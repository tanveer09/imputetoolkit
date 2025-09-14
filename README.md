# imputetoolkit

*Evaluate and compare multiple imputation methods with **consistent metrics** and an **intuitive S3 interface**.*

`imputetoolkit` is an R package for evaluating the quality of data imputation methods.\
This package implements an **object-oriented S3 interface** around an `evaluator` class that computes multiple metrics comparing imputed values with ground truth data.\
It is designed to help researchers and practitioners benchmark different imputation strategies side-by-side, providing both per-column and global metrics such as RMSE, MAE, R², correlation recovery, KS statistic, and accuracy.\
By wrapping results in an `evaluator` object, the package offers a consistent, user-friendly interface with familiar methods like `print()` and `summary()`.

------------------------------------------------------------------------

## Justification for OO Programming

The **evaluator** function is a strong candidate for Object-Oriented programming because:

1.  **Encapsulation of Related Components**
    -   Imputation evaluation naturally produces multiple outputs: per-variable metrics (e.g., RMSE, MAE, R²), global summaries, and metadata (method used, variables evaluated).\
    -   OO design allows these pieces to live inside a single structured object (`class = "evaluator"`) instead of scattering them across lists or separate return values.
2.  **Clear and Familiar Workflow for Users**
    -   Users can call familiar generic functions such as `print()`, `summary()`, or even `plot()` (if extended) on the evaluator object.\
    -   This design mirrors R’s built-in modeling ecosystem (`lm`, `glm`, `kmeans`), making it intuitive and lowering the learning curve.
3.  **Extensibility for Future Development**
    -   New metrics (e.g., correlation recovery, KS statistic) can be easily added without changing how users interact with the object.\
    -   Additional methods (`plot.evaluator`, `predict.evaluator`) could later extend functionality without rewriting core code.
4.  **Consistency Across Analyses**
    -   OO programming enforces a consistent structure: no matter which imputation method is evaluated, the output object behaves the same way.\
    -   This helps ensure reproducibility and comparability across datasets.
5.  **Alternatives Considered**
    -   A *functional* approach (returning a plain list) would work, but it would force the user to manually extract fields like `$rmse` or `$r2`, making the workflow clunkier.\
    -   R6 or Reference Classes could also be used, but for statistical models in R, **S3 classes are lightweight, idiomatic, and align with existing practices**.

Overall, the evaluator is a **good candidate for OO programming** because it bundles rich, structured outputs into an intuitive object, provides a user-friendly interface, and remains extensible for future enhancements.

------------------------------------------------------------------------

## Objects and Methods

-   **`evaluator()`**\
    Constructor that creates an object of class `"evaluator"`.\
    Takes two named lists of numeric vectors (`true_data`, `imputed_data`) and a method name.

-   **S3 Methods**

    -   `print.evaluator(x)` – displays global evaluation metrics for an imputation method.
    -   `summary.evaluator(x)` – returns a `data.frame` with per-column metrics and global averages.

-   **Metrics Computed** For each column:

    -   Root Mean Squared Error (RMSE)
    -   Mean Absolute Error (MAE)
    -   R² (coefficient of determination)
    -   Correlation recovery
    -   Kolmogorov–Smirnov statistic (distribution similarity)
    -   Accuracy (proportion of masked values correctly recovered)

    These are also aggregated into **global values** stored in the `evaluator` object.

------------------------------------------------------------------------

## Installation

You can install directly from GitHub:

``` r
# Install from GitHub
devtools::install_github("tanveer09/imputetoolkit")
```

------------------------------------------------------------------------

## Usage Example

``` r
library(imputetoolkit)

# Ground truth and imputed data
true_data <- list(
  age = c(25, 30, 40),
  income = c(50000, 60000, 70000)
)

imputed_data <- list(
  age = c(25, 31, 39),
  income = c(50000, 61000, 69000)
)

# Create evaluator object
result <- evaluator(true_data, imputed_data, method = "mean")

# Inspect results
print(result)
summary(result)
```

## Output (example)

```         
Evaluation for method: mean
Global Metrics:
  RMSE       : 1.2909
  MAE        : 1.0000
  R^2        : 0.9456
  Correlation: 0.9827
  KS         : 0.2000
  Accuracy   : 0.5000

Per-column metrics available in result$metrics
```

------------------------------------------------------------------------

## Workflow

1.  **Perform** multiple imputations on your dataset (mean, median, kNN, MICE, etc.).
2.  **Call** `evaluator(true_data, imputed_data, method)` for each method.
3.  **Collect** results into a list and compare across methods.
4.  **Use** `print()` for quick checks and `summary()` for detailed per-column analysis.

------------------------------------------------------------------------

## Testing

Unit tests are provided under the `tests/testthat/` directory. To run all tests:

``` r
devtools::test()
```

These tests check that:

-   Objects of class `"evaluator"` are created correctly.
-   Metrics are computed accurately for numeric data.
-   Errors are raised for invalid inputs (e.g., mismatched keys, NA/Inf values).
-   S3 methods (print, summary) return expected outputs.

------------------------------------------------------------------------

## Documentation

All functions are documented with **Roxygen2**. To rebuild documentation, run:

``` r
devtools::document()
```

Help pages are available for all major functions:

``` r
?evaluator
?print.evaluator
?summary.evaluator
```

------------------------------------------------------------------------

## LLM Disclosure

Some parts of this package — including **documentation drafting, README preparation, and sections of the R/C++ code (e.g., error handling and function scaffolding)**, were assisted by **ChatGPT (OpenAI)**.

All generated content was **reviewed, debugged, and adapted** before inclusion in the final submission.


