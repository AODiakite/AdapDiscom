# A multimodal high dimensionnal method accounting for missing data and
measurement error heterogenity


-   [AdapDiscom Package](#adapdiscom-package)
    -   [Installation](#installation)
    -   [Main Functions](#main-functions)
    -   [Basic Usage Example](#basic-usage-example)
        -   [Data Simulation](#data-simulation)
        -   [Running AdapDiscom](#running-adapdiscom)
        -   [Examining Results](#examining-results)
        -   [Comparing Methods](#comparing-methods)
    -   [Robust Estimation](#robust-estimation)
    -   [Different Covariance
        Structures](#different-covariance-structures)
    -   [Advanced Parameters](#advanced-parameters)
    -   [Block Structure Requirements](#block-structure-requirements)
    -   [Tips for Usage](#tips-for-usage)
    -   [References](#references)

# AdapDiscom Package

The `AdapDiscom` package provides implementations of a multimodal high
dimensionnal method accounting for missing data and measurement error
heterogenity

## Installation

``` r
# Install from local source
devtools::install_local("https://github.com/AODiakite/AdapDiscom")

# Load the package
library(AdapDiscom)
```

## Main Functions

The package provides four main estimation functions:

1.  **`adapdiscom()`**
2.  **`discom()`**  
3.  **`fast_adapdiscom()`**
4.  **`fast_discom()`**

## Basic Usage Example

### Data Simulation

Let’s start by simulating data with a block-diagonal covariance
structure:

``` r
library(AdapDiscom)
library(MASS)

# Set parameters
n <- 100        # Training sample size
n.tuning <- 50  # Tuning sample size  
n.test <- 50    # Test sample size
p <- 50         # Number of variables
pp <- c(20, 15, 15)  # Block sizes (must sum to p)

# Generate block-diagonal covariance matrix
Sigma <- generate.cov(p, example = 2)  # Block diagonal structure

# Generate true beta coefficients
beta.true <- c(rep(2, 10), rep(0, 10), rep(1, 10), rep(0, 5), rep(-1, 10), rep(0, 5))

# Generate training data
X.train <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
y.train <- X.train %*% beta.true + rnorm(n)
n1=n2=n3=n4=n%/%4
p1 = pp[1]
p2 = pp[2]
p3 = pp[3]
X.train[(n1 + 1):(n1 + n2), (p1 + p2 + 1):(p1 + p2 + p3)] <- NA
X.train[(n1 + n2 + 1):(n1 + n2 + n3), (p1 + 1):(p1 + p2)] <- NA
X.train[(n1 + n2 + n3 + 1):(n1 + n2 + n3 + n4), (1:p1)] <- NA

# Generate tuning data
X.tuning <- mvrnorm(n.tuning, mu = rep(0, p), Sigma = Sigma)
y.tuning <- X.tuning %*% beta.true + rnorm(n.tuning)

# Generate test data
X.test <- mvrnorm(n.test, mu = rep(0, p), Sigma = Sigma)
y.test <- X.test %*% beta.true + rnorm(n.test)
```

### Running AdapDiscom

``` r
# Run AdapDiscom with default parameters
result <- adapdiscom(
  beta = beta.true,      # True coefficients (optional, for evaluation)
  x = X.train,          # Training data
  y = y.train,          # Training response
  x.tuning = X.tuning,  # Tuning data
  y.tuning = y.tuning,  # Tuning response
  x.test = X.test,      # Test data
  y.test = y.test,      # Test response
  nlambda = 20,         # Number of lambda values
  nalpha = 10,          # Number of alpha values
  pp = pp,              # Block sizes
  robust = 0,           # Classical estimation
  standardize = TRUE,   # Standardize data
  itcp = TRUE          # Include intercept
)

# View results
print(result$R2)
```

              [,1]
    [1,] 0.8356919

### Examining Results

The functions return a list with the following components:

``` r
# Optimal parameters
cat("Optimal lambda:", result$lambda, "\n")
```

    Optimal lambda: 0.5802845 

``` r
cat("Optimal alpha:", result$alpha, "\n")
```

    Optimal alpha: 0.4641589 1 0.4641589 0.4641589 

``` r
# Performance metrics
cat("Training error:", result$train.error, "\n")
```

    Training error: 16.38544 

``` r
cat("Test error:", result$test.error, "\n")
```

    Test error: 15.15845 

``` r
cat("R-squared:", result$R2, "\n")
```

    R-squared: 0.8356919 

``` r
# Model selection
cat("Number of selected variables:", result$select, "\n")
```

    Number of selected variables: 35 

``` r
# If true beta provided, evaluation metrics
if (!is.null(beta.true)) {
  cat("Estimation error:", result$est.error, "\n")
  cat("False Positive Rate:", result$fpr, "\n")
  cat("False Negative Rate:", result$fnr, "\n")
}
```

    Estimation error: 3.722019 
    False Positive Rate: 0.35 
    False Negative Rate: 0.06666667 

``` r
# Estimated coefficients
head(result$a1)  # First 6 coefficients
```

    [1] 2.801564 1.836431 1.206341 1.242693 2.463127 1.915449

### Comparing Methods

``` r
# Compare different methods
methods <- c("adapdiscom", "discom", "fast_adapdiscom", "fast_discom")
results <- list()

# AdapDiscom
results$adapdiscom <- adapdiscom(
  beta = beta.true, x = X.train, y = y.train,
  x.tuning = X.tuning, y.tuning = y.tuning,
  x.test = X.test, y.test = y.test,
  nlambda = 20, nalpha = 10, pp = pp
)

# DISCOM
results$discom <- discom(
  beta = beta.true, x = X.train, y = y.train,
  x.tuning = X.tuning, y.tuning = y.tuning,
  x.test = X.test, y.test = y.test,
  nlambda = 20, nalpha = 10, pp = pp
)

# Fast AdapDiscom
results$fast_adapdiscom <- fast_adapdiscom(
  beta = beta.true, x = X.train, y = y.train,
  x.tuning = X.tuning, y.tuning = y.tuning,
  x.test = X.test, y.test = y.test,
  nlambda = 20, nalpha = 10, pp = pp
)

# Fast DISCOM
results$fast_discom <- fast_discom(
  beta = beta.true, x = X.train, y = y.train,
  x.tuning = X.tuning, y.tuning = y.tuning,
  x.test = X.test, y.test = y.test,
  nlambda = 20, pp = pp
)

# Compare performance
comparison <- data.frame(
  Method = names(results),
  Test_Error = sapply(results, function(x) x$test.error),
  R_Squared = sapply(results, function(x) x$R2),
  Selected_Vars = sapply(results, function(x) x$select),
  Est_Error = sapply(results, function(x) x$est.error),
  FPR = sapply(results, function(x) x$fpr),
  FNR = sapply(results, function(x) x$fnr),
  Time = sapply(results, function(x) x$time)
)

print(comparison)
```

                             Method Test_Error R_Squared Selected_Vars Est_Error
    adapdiscom           adapdiscom   15.15845 0.8356919            35  3.722019
    discom                   discom   10.41532 0.8858421            37  3.410328
    fast_adapdiscom fast_adapdiscom   12.97211 0.8387977            35  3.634373
    fast_discom         fast_discom   12.97211 0.8387977            35  3.634373
                     FPR        FNR   Time
    adapdiscom      0.35 0.06666667 54.408
    discom          0.40 0.03333333  0.386
    fast_adapdiscom 0.40 0.10000000  0.244
    fast_discom     0.40 0.10000000  0.082

## Robust Estimation

The package supports robust estimation using Huber’s M-estimator:

``` r
# Run with robust estimation
result_robust <- adapdiscom(
  beta = beta.true,
  x = X.train, y = y.train,
  x.tuning = X.tuning, y.tuning = y.tuning,
  x.test = X.test, y.test = y.test,
  nlambda = 20, nalpha = 10, pp = pp,
  robust = 1,        # Enable robust estimation
  k.value = 1.5      # Huber tuning parameter
)
```

## Different Covariance Structures

The package handles different covariance structures:

``` r
# AR(1) structure
Sigma_ar1 <- generate.cov(p, example = 1)

# Block diagonal structure
Sigma_block <- generate.cov(p, example = 2)

# Kronecker product structure
Sigma_kron <- generate.cov(p, example = 3)
```

## Advanced Parameters

``` r
# Advanced usage with custom parameters
result_advanced <- adapdiscom(
  beta = beta.true,
  x = X.train, y = y.train,
  x.tuning = X.tuning, y.tuning = y.tuning,
  x.test = X.test, y.test = y.test,
  nlambda = 30,              # More lambda values
  nalpha = 15,               # More alpha values
  pp = pp,
  robust = 0,
  standardize = TRUE,
  itcp = TRUE,
  lambda.min.ratio = 0.001,  # Custom lambda range
  k.value = 1.5
)
```

## Block Structure Requirements

-   **AdapDiscom**: Supports any number of blocks (specified by `pp`
    vector)
-   **DISCOM**: Supports 2, 3, or 4 blocks only
-   **Fast variants**: Optimized versions with reduced computational
    complexity

## Tips for Usage

1.  **Block sizes**: Ensure `sum(pp) == p` (total number of variables)
2.  **Sample sizes**: Use adequate tuning and test sample sizes for
    reliable parameter selection
3.  **Standardization**: Generally recommended (`standardize = TRUE`)
4.  **Robust estimation**: Use when data contains outliers
    (`robust = 1`)
5.  **Lambda range**: Adjust `lambda.min.ratio` if needed for better
    regularization path

## References

For more details on the methodology, please refer to the original papers
on *A multimodal high dimensionnal method accounting for missing data
and measurement error heterogenity*.
