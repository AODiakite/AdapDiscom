#' Generate Covariance Matrix
#'
#' @param p Integer, dimension of the covariance matrix
#' @param example Integer, type of covariance structure (1=AR(1), 2=Block diagonal, 3=Kronecker product)
#' @return A p x p covariance matrix
#' @importFrom stats cor cov na.omit
#' @examples
#' # AR(1) covariance structure
#' Sigma1 <- generate.cov(p = 20, example = 1)
#' print(Sigma1[1:3, 1:3])
#' 
#' # Block diagonal structure (p must be multiple of 5)
#' Sigma2 <- generate.cov(p = 25, example = 2) 
#' print(Sigma2[1:5, 1:5])
#' 
#' # Kronecker product structure (p must be multiple of 10)
#' Sigma3 <- generate.cov(p = 100, example = 3)
#' print(dim(Sigma3))
#' @export
generate.cov <- function(p, example) {
  cov.matrix <- matrix(0, p, p)
  if (example == 1) {
    # AR(1) structure
    for (i in 1:p) {
      for (j in 1:p) {
        cov.matrix[i, j] <- 0.6^(abs(i - j))
      }
    }
  } else if (example == 2) {
    # Block diagonal structure
    for (i in 1:(p/5)) {
      cov.matrix[(1 + (i-1)*5):(5*i), (1 + (i-1)*5):(5*i)] <- matrix(0.15, 5, 5) + 0.85*diag(5)
    }
  } else {
    # Pour cov.matrix1 (10×10)
    cov.matrix1 <- matrix(NA, 10, 10)
    # Remplir avec des valeurs de corrélation appropriées
    for (i in 1:10) {
      for (j in 1:10) {
        if (i == j) cov.matrix1[i, j] <- 1
        else cov.matrix1[i, j] <- 0.8^(abs(i - j))
      }
    }
    
    # Pour cov.matrix2 (10×10)
    pnew <- p/10
    cov.matrix2 <- matrix(NA, pnew, pnew)
    for (i in 1:pnew) {
      for (j in 1:pnew) {
        cov.matrix2[i, j] <- 0.3^(abs(i - j))
      }
    }
    
    # Produit de Kronecker
    cov.matrix <- kronecker(cov.matrix1, cov.matrix2)
  }
  cov.matrix
}

#' Compute X'X Matrix
#'
#' @param x Matrix, input data matrix
#' @param robust Integer, 0 for classical estimate, 1 for Huber robust estimate
#' @param k_value Numeric, tuning parameter for Huber function
#' @return Covariance matrix
#' @examples
#' # Create sample data with missing values
#' set.seed(123)
#' x <- matrix(rnorm(100), 20, 5)
#' x[sample(100, 10)] <- NA  # Introduce missing values
#' 
#' # Classical covariance estimation
#' xtx_classical <- compute.xtx(x, robust = 0)
#' print(round(xtx_classical, 3))
#' 
#' # Robust covariance estimation
#' xtx_robust <- compute.xtx(x, robust = 1, k_value = 1.5)
#' print(round(xtx_robust, 3))
#' @export
compute.xtx <- function(x, robust = 0, k_value = 1.5) {
  p <- ncol(x)
  cov.matrix <- matrix(NA, p, p)
  if (robust == 0) {
    cov.matrix <- cov(x, use = "pairwise.complete.obs")
  } else {
    for (i in 1:p) {
      for (j in 1:p) {
        index <- which((!is.na(x[, i])) & (!is.na(x[, j])))
        x1 <- x[index, i]
        x2 <- x[index, j]
        cov.matrix[i, j] <- robustbase::huberM((x1 - mean(x1)) * (x2 - mean(x2)), 
                                               k = k_value * sqrt(length(index)/log(p)))$mu
      }
    }
  }
  cov.matrix
}

#' Compute X'y Vector
#'
#' @param x Matrix, input data matrix
#' @param y Vector, response vector
#' @param robust Integer, 0 for classical estimate, 1 for Huber robust estimate
#' @param k_value Numeric, tuning parameter for Huber function
#' @return Covariance vector
#' @examples
#' # Create sample data
#' set.seed(123)
#' x <- matrix(rnorm(100), 20, 5)
#' y <- rnorm(20)
#' x[sample(100, 8)] <- NA  # Missing values in x
#' 
#' # Classical cross-covariance
#' xty_classical <- compute.xty(x, y, robust = 0)
#' print(round(xty_classical, 3))
#' 
#' # Robust cross-covariance
#' xty_robust <- compute.xty(x, y, robust = 1)
#' print(round(xty_robust, 3))
#' @export
compute.xty <- function(x, y, robust = 0, k_value = 1.5) {
  p <- ncol(x)
  cov.vector <- rep(NA, p)
  if (robust == 0) {
    cov.vector <- cov(x, y, use = "pairwise.complete.obs")
  } else {
    for (i in 1:p) {
      index <- which((!is.na(x[, i])))
      x1 <- x[index, i]
      x2 <- y[index]
      cov.vector[i] <- robustbase::huberM((x1 - mean(x1)) * (x2 - mean(x2)), 
                                          k = k_value * sqrt(length(index)/log(p)))$mu
    }
  }
  cov.vector
}

standardize_x <- function(x, robust = 0, k.value = 1.5) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  x.mean <- colMeans(x, na.rm = TRUE)
  x.sd <- sqrt(diag(compute.xtx(x, robust = robust, k.value)))
  x.sd[x.sd < 1e-6] <- 1
  x <- sweep(x, 2, x.mean, '-')
  x <- sweep(x, 2, x.sd, '/')
  x <- as.matrix(x)
  
  return(list(x = x, x.mean = x.mean, x.sd = x.sd))
}

fit_standardize_x <- function(x, x.mean, x.sd) {
  x <- as.matrix(x)
  x <- sweep(x, 2, x.mean, '-')
  x <- sweep(x, 2, x.sd, '/')
  return(as.matrix(x))
}

#' Compute Lambda Max for L1 Regularization using KKT Conditions
#'
#' @param X Matrix, design matrix
#' @param y Vector, response vector
#' @param Methode Character, method for computation
#' @param robust Integer, 0 for classical, 1 for robust
#' @return Maximum lambda value
#' @examples
#' # Generate sample data
#' set.seed(123)
#' n <- 50; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' y <- rnorm(n)
#' 
#' # Different methods for lambda_max computation
#' lmax_lasso <- lambda_max(X, y, Methode = "lasso")
#' lmax_discom <- lambda_max(X, y, Methode = "discom")
#' 
#' print(paste("Lambda max (lasso):", round(lmax_lasso, 4)))
#' print(paste("Lambda max (discom):", round(lmax_discom, 4)))
#' @export
lambda_max <- function(X, y, Methode = "lasso", robust = 0) {
  n <- nrow(X)
  if (Methode == "lasso") {
    data <- cbind(y, X)
    data <- na.omit(data)
    X <- as.matrix(data[, -1])
    y <- as.matrix(data[, 1])
    lambda_max <- t(X) %*% y
    lambda_max <- max(abs(lambda_max)) / n
  } else if (Methode == "mean") {
    X <- im_mean(X)
    lambda_max <- t(X) %*% y
    lambda_max <- max(abs(lambda_max)) / n
  } else if (Methode == "svd") {
    X <- im_svd(X)
    lambda_max <- t(X) %*% y
    lambda_max <- max(abs(lambda_max)) / n
  } else {
    Xy <- compute.xty(X, y, robust = robust)
    lambda_max <- max(abs(Xy))
  }
  return(lambda_max)
}

#' Get Block Indices
#'
#' @param pp Vector, block sizes
#' @return List with start and end indices for each block
#' @examples
#' # Define block sizes
#' pp <- c(10, 15, 20)
#' indices <- get_block_indices(pp)
#' print(indices)
#' # Shows: $starts = c(1, 11, 26) and $ends = c(10, 25, 45)
#' 
#' # For two blocks
#' pp2 <- c(25, 25)
#' indices2 <- get_block_indices(pp2)
#' print(indices2)
#' @export
get_block_indices <- function(pp) {
  n <- length(pp)
  starts <- c(1, cumsum(pp[-n]) + 1)
  ends <- cumsum(pp)
  return(list(starts = starts, ends = ends))
}

im_mean <- function(X) {
  X <- as.matrix(X)
  for (j in 1:ncol(X)) {
    missing_idx <- is.na(X[, j])
    if (any(missing_idx)) {
      X[missing_idx, j] <- mean(X[, j], na.rm = TRUE)
    }
  }
  return(X)
}

im_svd <- function(X) {
  X <- as.matrix(X)
  if (any(is.na(X))) {
    # Use softImpute for SVD-based imputation
    result <- softImpute::softImpute(X, lambda = 0, maxit = 100)
    X <- result$u %*% diag(result$d) %*% t(result$v)
  }
  return(X)
}

