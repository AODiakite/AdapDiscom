# A multimodal high dimensionnal method accounting for missing data and
measurement error heterogenity


-   [AdapDiscom Package](#adapdiscom-package)
    -   [Installation](#installation)
    -   [Main Functions](#main-functions)
    -   [Different Covariance
        Structures](#different-covariance-structures)
    -   [Basic Usage Example](#basic-usage-example)
        -   [Data Simulation](#data-simulation)
        -   [Running AdapDiscom](#running-adapdiscom)
        -   [Examining Results](#examining-results)
        -   [Comparing Methods](#comparing-methods)
    -   [Robust Estimation](#robust-estimation)
    -   [Advanced Parameters](#advanced-parameters)
    -   [Block Structure Requirements](#block-structure-requirements)
    -   [Tips for Usage](#tips-for-usage)
    -   [References](#references)

# AdapDiscom Package

The `AdapDiscom` package provides implementations of **AdapDISCOM**, a
novel adaptive direct sparse regression method for high-dimensional
multimodal data with block-wise missingness and measurement errors.
Building on the DISCOM framework, AdapDISCOM introduces
modality-specific weighting schemes to account for heterogeneity in data
structures and error magnitudes across modalities. Unlike existing
methods that are limited to a fixed number of blocks, **AdapDISCOM
supports any number of blocks K**, making it highly flexible for diverse
multimodal data applications. The package includes robust variants
(AdapDISCOM-Huber) and computationally efficient versions
(Fast-AdapDISCOM) to handle heavy-tailed distributions and large-scale
datasets. AdapDISCOM has been shown to consistently outperform existing
methods such as DISCOM, SCOM, and imputation based, particularly under
heterogeneous contamination scenarios, providing a scalable framework
for realistic multimodal data analysis with missing data and measurement
errors.

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

## Different Covariance Structures

The package handles different covariance structures:

``` r
# Number of variables
p <-  10
# Load the package
library(AdapDiscom)
# AR(1) structure
Sigma_ar1 <- generate.cov(p, example = 1)

# Block diagonal structure
Sigma_block <- generate.cov(p, example = 2)

# Kronecker product structure
Sigma_kron <- generate.cov(p, example = 3)
```

## Basic Usage Example

### Data Simulation

This simulation scenario represents a simple replication of scenario 2
from the AdapDISCOM’s paper. It is case study with 300 variables
distributed across 3 equal blocks (100 variables each), where each block
follows a different covariance structure: the first block uses an AR(1)
structure, the second a block-diagonal structure, and the third a
Kronecker product structure. The training data (n=440) exhibits a
structured missing data pattern where each quarter of the sample is
completely missing one of the three variable blocks, creating
heterogeneity in data availability that reflects real-world situations
where different subgroups of subjects may have access to different sets
of measurements. The true coefficient vector *β* follows a sparse
pattern with only 5% of variables having non-zero effects (*β* = 0.5),
allowing evaluation of the method’s ability to correctly identify
relevant variables in a high-dimensional context with heterogeneous
missing data.

``` r
library(AdapDiscom)
library(MASS)
set.seed(123)
# Set parameters
p <- 300
n <- 440
n.tuning <- 200
n.test <- 400
p1 <- p%/%3
p2 <- p%/%3
p3 <- p - p1 - p2
pp <- c(p1, p2, p3)  # Block sizes
sigma <- 1
```

``` r
# Generate different covariance matrices for each block
cov.mat1 <- generate.cov(p1, 1)  # AR(1) structure for block 1
cov.mat2 <- generate.cov(p2, 2)  # Block diagonal for block 2
cov.mat3 <- generate.cov(p3, 3)  # Kronecker product for block 3
```

``` r
# Generate true beta coefficients
beta1 <- rep(c(rep(0.5, 5), rep(0, 95)), p/100)
beta.true <- beta1

# Generate complete data for all samples
pre.x1 <- mvrnorm(n + n.tuning + n.test, rep(0, p%/%3), cov.mat1)
pre.x2 <- mvrnorm(n + n.tuning + n.test, rep(0, p%/%3), cov.mat2)
pre.x3 <- mvrnorm(n + n.tuning + n.test, rep(0, p%/%3), cov.mat3)
pre.x <- cbind(pre.x1, pre.x2, pre.x3)
pre.ep <- rnorm(n + n.tuning + n.test, 0, sigma)

# Training data
n.com <- n/4
n1 <- n2 <- n3 <- n4 <- n.com
X.train <- pre.x[1:n, ]
ep <- pre.ep[1:n]
y.train <- X.train %*% beta1 + ep 
colnames(X.train) <- paste0("X", 1:p)

# Introduce missing data pattern
X.train[(n1 + 1):(n1 + n2), (p1 + p2 + 1):(p1 + p2 + p3)] <- NA
X.train[(n1 + n2 + 1):(n1 + n2 + n3), (p1 + 1):(p1 + p2)] <- NA
X.train[(n1 + n2 + n3 + 1):(n1 + n2 + n3 + n4), (1:p1)] <- NA

# Tuning data
X.tuning <- pre.x[(n + 1):(n + n.tuning), ]
ep.tuning <- pre.ep[(n + 1):(n + n.tuning)]
y.tuning <- X.tuning %*% beta1 + ep.tuning
colnames(X.tuning) <- paste0("X", 1:p)

# Test data
X.test <- pre.x[(n + n.tuning + 1):(n + n.tuning + n.test), ]
ep.test <- pre.ep[(n + n.tuning + 1):(n + n.tuning + n.test)]
y.test <- X.test %*% beta1 + ep.test
colnames(X.test) <- paste0("X", 1:p)
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
    [1,] 0.8089201

### Examining Results

The functions return a list with the following components:

``` r
# Optimal parameters
cat("Optimal lambda:", result$lambda, "\n")
```

    Optimal lambda: 0.2823072 

``` r
cat("Optimal alpha:", result$alpha, "\n")
```

    Optimal alpha: 0.4641589 0.2154435 0.4641589 0.4641589 

``` r
# Performance metrics
cat("Training error:", result$train.error, "\n")
```

    Training error: 1.178665 

``` r
cat("Test error:", result$test.error, "\n")
```

    Test error: 1.526394 

``` r
cat("R-squared:", result$R2, "\n")
```

    R-squared: 0.8089201 

``` r
# Model selection
cat("Number of selected variables:", result$select, "\n")
```

    Number of selected variables: 27 

``` r
# If true beta provided, evaluation metrics
if (!is.null(beta.true)) {
  cat("Estimation error:", result$est.error, "\n")
  cat("False Positive Rate:", result$fpr, "\n")
  cat("False Negative Rate:", result$fnr, "\n")
}
```

    Estimation error: 0.7911645 
    False Positive Rate: 0.04210526 
    False Negative Rate: 0 

``` r
# Estimated coefficients
head(result$a1)  # First 6 coefficients
```

    [1] 0.36421275 0.57343638 0.70592242 0.56074589 0.42810130 0.09403387

### Comparing Methods

``` r
# Compare different methods
methods <- c("adapdiscom", "discom", "fast_adapdiscom", "fast_discom")
results <- list()

# AdapDiscom
results$adapdiscom <- result # Already computed above

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
    adapdiscom           adapdiscom   1.526394 0.8089201            27 0.7911645
    discom                   discom   1.698120 0.7851645            27 0.8079454
    fast_adapdiscom fast_adapdiscom   1.471665 0.8173831            33 0.7132907
    fast_discom         fast_discom   1.471665 0.8173831            33 0.7132907
                           FPR FNR    Time
    adapdiscom      0.04210526   0 465.022
    discom          0.04210526   0   5.422
    fast_adapdiscom 0.06315789   0   1.825
    fast_discom     0.06315789   0   0.897

> **Note**
>
> When the missing blocks are of equal size within each modality and are
> disjoint across modalities (i.e., no observations have missing data in
> multiple modalities simultaneously), Fast-AdapDISCOM reduces to
> Fast-DISCOM since the adaptive weighting scheme becomes uniform across
> all modalities.

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
