#' DISCOM: Optimal Sparse Linear Prediction for Block-missing Multi-modality Data Without Imputation
#'
#' @param beta Vector, true beta coefficients (optional)
#' @param x Matrix, training data
#' @param y Vector, training response
#' @param x.tuning Matrix, tuning data
#' @param y.tuning Vector, tuning response
#' @param x.test Matrix, test data
#' @param y.test Vector, test response
#' @param nlambda Integer, number of lambda values
#' @param nalpha Integer, number of alpha values
#' @param pp Vector, block sizes. Discom supports 2, 3, or 4 blocks.
#' @param robust Integer, 0 for classical, 1 for robust estimation
#' @param standardize Logical, whether to standardize covariates. When TRUE, uses training data mean and standard deviation to standardize tuning and test sets. When robust=1, uses Huber-robust standard deviation estimates
#' @param itcp Logical, whether to include intercept
#' @param lambda.min.ratio Numeric, `lambda.min.ratio` sets the smallest lambda value in the grid, expressed as a fraction of `lambda.max`—the smallest lambda for which all coefficients are zero. By default, it is `0.0001` when the number of observations (`nobs`) exceeds the number of variables (`nvars`), and `0.01` when `nobs < nvars`. Using a very small value in the latter case can lead to overfitting.
#' @param k.value Numeric, tuning parameter for robust estimation
#'
#' @return List with estimation results
#' @section Value:
#' The function returns a list containing the following components:
#'
#' \describe{
#'   \item{err}{A multi-dimensional array storing the mean squared error (MSE) for all combinations of tuning parameters alpha and lambda.}
#'   \item{est.error}{The estimation error, calculated as the Euclidean distance between the estimated beta coefficients and the true beta (if provided).}
#'   \item{lambda}{The optimal lambda value chosen via cross-validation on the tuning set.}
#'   \item{alpha}{A vector of the optimal alpha values, also selected on the tuning set.}
#'   \item{train.error}{The mean squared error on the tuning set for the optimal parameter combination.}
#'   \item{test.error}{The mean squared error on the test set for the final, optimal model.}
#'   \item{y.pred}{The predicted values for the observations in the test set.}
#'   \item{R2}{The R-squared value, which measures the proportion of variance explained by the model on the test set.}
#'   \item{a0}{The intercept of the final model.}
#'   \item{a1}{The vector of estimated beta coefficients for the final model.}
#'   \item{select}{The number of non-zero coefficients, representing the number of selected variables.}
#'   \item{xtx}{The final regularized covariance matrix used to fit the optimal model.}
#'   \item{fpr}{The False Positive Rate (FPR) if the true beta is provided. It measures the proportion of irrelevant variables incorrectly selected.}
#'   \item{fnr}{The False Negative Rate (FNR) if the true beta is provided. It measures the proportion of relevant variables incorrectly excluded.}
#'   \item{lambda.all}{The complete vector of all lambda values tested during cross-validation.}
#'   \item{beta.cov.lambda.max}{The estimated beta coefficients using the maximum lambda value.}
#'   \item{time}{The total execution time of the function in seconds.}
#' }
#' @export
discom <- function(beta, x, y, x.tuning, y.tuning, x.test, y.test, nlambda, nalpha, pp, 
                   robust = 0, standardize = TRUE, itcp = TRUE, lambda.min.ratio = NULL,
                   k.value = 1.5) {
  
  x_std = standardize_x(as.matrix(x), robust = robust, k.value = k.value)
  xm <- x_std$x.mean
  xs <- x_std$x.sd
  if (standardize) {
    x <- x_std$x
    x.tuning <- fit_standardize_x(x.tuning, xm, xs)
    x.test <- fit_standardize_x(x.test, xm, xs)
  }
  
  # Temps de calcul
  start_time <- proc.time()
  n <- dim(x)[1]
  p <- dim(x)[2]
  n.tuning <- dim(x.tuning)[1]
  n.test <- dim(x.test)[1]
  
  alpha.all <- 10^seq(0, -3, length = nalpha)
  lmax <- lambda_max(x, y, Methode = "discom", robust = robust)
  nobs <- dim(na.omit(x))[1]
  if(is.null(lambda.min.ratio) || length(lambda.min.ratio) == 0 || lambda.min.ratio <= 0 ||
     lambda.min.ratio > 1){
    lambda.min.ratio <- ifelse(nobs < p, 0.01, 1e-04)
  }
  
  lmin <- lmax * lambda.min.ratio
  lambda.all <- exp(seq(log(lmax), log(lmin), length.out = nlambda))
  
  DISCOM.tun.error <- array(NA, dim = c(nalpha, nalpha, nlambda))
  
  xtx.raw <- compute.xtx(x, robust = robust, k_value = k.value)
  xty <- compute.xty(x, y, robust = robust, k_value = k.value)
  
  if (length(pp) == 4) {
    p1 <- pp[1]
    p2 <- pp[2]
    p3 <- pp[3]
    p4 <- pp[4]
    xtx.raw.I <- as.matrix(Matrix::bdiag(xtx.raw[1:p1, 1:p1], 
                                        xtx.raw[(p1+1):(p1+p2), (p1+1):(p1+p2)],
                                        xtx.raw[(p1+p2+1):(p1+p2+p3), (p1+p2+1):(p1+p2+p3)],
                                        xtx.raw[(p1+p2+p3+1):(p1+p2+p3+p4), (p1+p2+p3+1):(p1+p2+p3+p4)]))
    xtx.raw.C <- xtx.raw - xtx.raw.I
    shrink.target <- sum(diag(xtx.raw))/p
  } else if (length(pp) == 3) {
    p1 <- pp[1]
    p2 <- pp[2]
    p3 <- pp[3]
    xtx.raw.I <- as.matrix(Matrix::bdiag(xtx.raw[1:p1, 1:p1], 
                                        xtx.raw[(p1+1):(p1+p2), (p1+1):(p1+p2)],
                                        xtx.raw[(p1+p2+1):(p1+p2+p3), (p1+p2+1):(p1+p2+p3)]))
    xtx.raw.C <- xtx.raw - xtx.raw.I
    shrink.target <- sum(diag(xtx.raw))/p
  } else if (length(pp) == 2) {
    p1 <- pp[1]
    p2 <- pp[2]
    xtx.raw.I <- as.matrix(Matrix::bdiag(xtx.raw[1:p1, 1:p1], 
                                        xtx.raw[(p1+1):(p1+p2), (p1+1):(p1+p2)]))
    xtx.raw.C <- xtx.raw - xtx.raw.I
    shrink.target <- sum(diag(xtx.raw))/p
  }
  
  for (i in 1:nalpha) {
    for (j in 1:nalpha) {
      alpha1 <- alpha.all[i]
      alpha2 <- alpha.all[j]
      xtx <- alpha1*xtx.raw.I + alpha2*xtx.raw.C + (1-alpha1)*shrink.target*diag(p)
      if (min(eigen(xtx)$values) < 0) {
        DISCOM.tun.error[i, j, ] <- 10^8
      } else {
        beta.initial <- rep(0, p)
        for (k in 1:nlambda) {
          beta.cov <- as.vector(scout::crossProdLasso(xtx, xty, lambda.all[k], 
                                                       beta.init = beta.initial)$beta)
          beta.initial <- beta.cov
          intercept <- ifelse(itcp, mean(y) - sum(beta.cov * (xm/xs)), 0)
          DISCOM.tun.values <- as.vector(as.matrix(x.tuning) %*% beta.cov) + intercept
          DISCOM.tun.error[i, j, k] <- mean((y.tuning - DISCOM.tun.values)^2)
        }
      }
    }
  }
  
  opt.index <- as.vector(which(DISCOM.tun.error == min(DISCOM.tun.error), arr.ind = TRUE)[1, ])
  train.error <- min(DISCOM.tun.error)
  opt.alpha1 <- alpha.all[opt.index[1]]
  opt.alpha2 <- alpha.all[opt.index[2]]
  opt.lambda <- lambda.all[opt.index[3]]
  xtx <- opt.alpha1*xtx.raw.I + opt.alpha2*xtx.raw.C + (1-opt.alpha1)*shrink.target*diag(p)
  beta.cov <- as.vector(scout::crossProdLasso(xtx, xty, opt.lambda)$beta)
  beta.cov.lambda.max <- as.vector(scout::crossProdLasso(xtx, xty, lmax)$beta)
  
  intercept <- ifelse(itcp, mean(y) - sum(beta.cov * (xm/xs)), 0)
  predict.test.values <- as.vector(x.test %*% beta.cov) + intercept
  
  DISCOM.test.error <- mean((y.test - predict.test.values)^2)
  select <- sum(as.vector(as.integer(beta.cov != 0)))
  
  if (!is.null(beta)) {
    DISCOM.fpr <- sum((beta == 0) & (beta.cov != 0))/sum(beta == 0)
    DISCOM.fnr <- sum((beta != 0) & (beta.cov == 0))/sum(beta != 0)
    DISCOM.est.error <- sqrt(sum((beta.cov - beta)^2))
  } else {
    DISCOM.fpr <- NA
    DISCOM.fnr <- NA
    DISCOM.est.error <- NA
  }
  
  R2 <- cor(predict.test.values, y.test)^2
  
  end_time <- proc.time()
  time_taken <- end_time - start_time
  
  a <- list('err' = DISCOM.tun.error, 'est.error' = DISCOM.est.error, 'lambda' = opt.lambda, 
           'alpha' = c(opt.alpha1, opt.alpha2), 'train.error' = train.error, 
           'test.error' = DISCOM.test.error, 'y.pred' = predict.test.values, 'R2' = R2,
           'lambda.all' = lambda.all, 'beta.cov.lambda.max' = beta.cov.lambda.max,
           'a0' = intercept, 'a1' = beta.cov, 'select' = select, 'xtx' = xtx, 
           'fpr' = DISCOM.fpr, 'fnr' = DISCOM.fnr, "time" = as.numeric(time_taken[3]))
  
  return(a)
}