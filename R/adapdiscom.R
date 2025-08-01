#' AdapDiscom: An Adaptive Sparse Regression Method for High-Dimensional Multimodal Data With Block-Wise Missingness and Measurement Errors
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
#' @param pp Vector, block sizes
#' @param robust Integer, 0 for classical, 1 for robust estimation of covariance
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
#' @examples
#' \dontrun{
#' # Generate synthetic data with block structure
#' set.seed(123)
#' n <- 100
#' pp <- c(20, 30)  # Two blocks
#' p <- sum(pp)
#' beta_true <- c(rep(0.5, 10), rep(0, 10), rep(0.3, 10), rep(0, 20))
#' 
#' # Generate covariance matrix
#' Sigma <- generate.cov(p = p, example = 1)
#' 
#' # Simulate data with missing blocks
#' x_full <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
#' # Introduce block-wise missing pattern
#' missing_pattern <- sample(1:2, n, replace = TRUE, prob = c(0.3, 0.7))
#' x <- x_full
#' x[missing_pattern == 1, (pp[1] + 1):p] <- NA  # Missing second block
#' 
#' # Generate response
#' y <- x_full %*% beta_true + rnorm(n, 0, 0.1)
#' 
#' # Split data
#' train_idx <- 1:60
#' tune_idx <- 61:80
#' test_idx <- 81:100
#' 
#' result <- adapdiscom(
#'   beta = beta_true,
#'   x = x[train_idx, ],
#'   y = y[train_idx],
#'   x.tuning = x[tune_idx, ],
#'   y.tuning = y[tune_idx],
#'   x.test = x[test_idx, ],
#'   y.test = y[test_idx],
#'   nlambda = 10,
#'   nalpha = 5,
#'   pp = pp
#' )
#' 
#' # View results
#' print(paste("Test error:", round(result$test.error, 4)))
#' print(paste("R-squared:", round(result$R2, 4)))
#' print(paste("Variables selected:", result$select))
#' }
#' @export
adapdiscom <- function(beta, x, y, x.tuning, y.tuning, x.test, y.test, nlambda, nalpha, pp, 
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
  
  start_time <- proc.time()
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  n.tuning <- dim(x.tuning)[1]
  n.test <- dim(x.test)[1]
  n_blocks <- length(pp)
  
  alpha.all <- 10^seq(10^(-10), -3, length = nalpha)
  
  lmax <- lambda_max(x, y, Methode = "discom", robust = robust)
  nobs <- dim(na.omit(x))[1]
  if(is.null(lambda.min.ratio) || length(lambda.min.ratio) == 0 || lambda.min.ratio <= 0 ||
     lambda.min.ratio > 1){
    lambda.min.ratio <- ifelse(nobs < p, 0.01, 1e-04)
  }
  
  lmin <- lmax * lambda.min.ratio
  lambda.all <- exp(seq(log(lmax), log(lmin), length.out = nlambda))
  
  # Création d'un tableau multidimensionnel pour les erreurs
  dim_array <- rep(nalpha, n_blocks + 1)
  dim_array <- c(dim_array, nlambda)
  DISCOM.tun.error <- array(NA, dim = dim_array)
  
  xtx.raw <- compute.xtx(x, robust = robust, k_value = k.value)
  xty <- compute.xty(x, y, robust = robust, k_value = k.value)
  
  # Création des indices pour les sous-matrices
  indices <- get_block_indices(pp)
  
  # Création des sous-matrices individuelles
  xtx.raw.blocks <- list()
  shrink.targets <- numeric(n_blocks)
  
  for (i in 1:n_blocks) {
    idx_range <- indices$starts[i]:indices$ends[i]
    xtx.raw.blocks[[i]] <- xtx.raw[idx_range, idx_range]
    shrink.targets[i] <- sum(diag(xtx.raw.blocks[[i]]))/p
  }
  
  # Création de la matrice bloc-diagonale initiale
  xtx.raw.I <- as.matrix(do.call(Matrix::bdiag, xtx.raw.blocks))
  
  # Calcul de xtx.raw.C
  xtx.raw.C <- xtx.raw - xtx.raw.I
  
  # Calcul du shrink.target global
  shrink.target <- sum(diag(xtx.raw))/p
  
  # Créer des indices pour les boucles imbriquées
  alpha_indices <- expand.grid(replicate(n_blocks + 1, 1:nalpha, simplify = FALSE))
  
  # Boucle sur toutes les combinaisons d'indices alpha
  for (row in 1:nrow(alpha_indices)) {
    # Extraire les valeurs alpha actuelles
    current_alphas <- alpha.all[as.numeric(alpha_indices[row, ])]
    
    # Les n_blocks premières alphas sont pour les blocs diagonaux
    block_alphas <- current_alphas[1:n_blocks]
    
    # Le dernier alpha est pour xtx.raw.C
    alpha_non_diag <- current_alphas[n_blocks + 1]
    
    # Calcul du shrink.target pondéré
    weighted_targets <- 0
    weight_sum <- 0
    for (i in 1:n_blocks) {
      weighted_targets <- weighted_targets + ((1-block_alphas[i])^2) * shrink.targets[i]
      weight_sum <- weight_sum + (1-block_alphas[i])^2
    }
    shrink.target_weighted <- weighted_targets / weight_sum
    
    # Construction de la matrice bloc-diagonale avec les alphas appliqués
    scaled_blocks <- list()
    for (i in 1:n_blocks) {
      scaled_blocks[[i]] <- block_alphas[i] * xtx.raw.blocks[[i]]
    }
    xtx.raw.I_scaled <- as.matrix(do.call(Matrix::bdiag, scaled_blocks))
    
    # Construction de la matrice finale
    xtx <- xtx.raw.I_scaled + alpha_non_diag * xtx.raw.C + 
           (n_blocks - sum(block_alphas)) * shrink.target_weighted * diag(p)
    
    # Vérification de la positivité
    if (min(eigen(xtx, only.values = TRUE)$values) < 0) {
      # Remplir avec 10^8 pour cette combinaison d'alphas et tous les lambdas
      idx <- as.list(as.numeric(alpha_indices[row, ]))
      idx_str <- paste(paste(idx, collapse = ","), ",", sep = "")
      eval(parse(text = paste("DISCOM.tun.error[", idx_str, "] = 10^8", sep = "")))
    } else {
      # Calcul des erreurs pour chaque lambda
      beta.initial <- rep(0, p)
      for (k in 1:nlambda) {
        beta.cov <- as.vector(scout::crossProdLasso(xtx, xty, lambda.all[k], 
                                                     beta.init = beta.initial)$beta)
        beta.initial <- beta.cov
        intercept <- ifelse(itcp, mean(y) - sum(beta.cov * (xm/xs)), 0)
        DISCOM.tun.values <- as.vector(as.matrix(x.tuning) %*% beta.cov) + intercept
        
        # Indexer correctement le tableau multidimensionnel
        idx <- c(as.numeric(alpha_indices[row, ]), k)
        idx_str <- paste(paste(idx, collapse = ","), sep = "")
        eval(parse(text = paste("DISCOM.tun.error[", idx_str, "] = mean((y.tuning-DISCOM.tun.values)^2)", sep = "")))
      }
    }
  }
  
  # Trouver les paramètres optimaux
  opt.index <- as.vector(which(DISCOM.tun.error == min(DISCOM.tun.error), arr.ind = TRUE)[1, ])
  train.error <- min(DISCOM.tun.error)
  
  # Extraire les valeurs alpha optimales
  opt.alphas <- alpha.all[opt.index[1:(n_blocks+1)]]
  block_alphas_opt <- opt.alphas[1:n_blocks]
  alpha_non_diag_opt <- opt.alphas[n_blocks+1]
  
  # Extraire le lambda optimal
  opt.lambda <- lambda.all[opt.index[n_blocks+2]]
  
  # Recalcul du shrink.target pondéré avec les alphas optimaux
  weighted_targets <- 0
  weight_sum <- 0
  for (i in 1:n_blocks) {
    weighted_targets <- weighted_targets + ((1-block_alphas_opt[i])^2) * shrink.targets[i]
    weight_sum <- weight_sum + (1-block_alphas_opt[i])^2
  }
  shrink.target_opt <- weighted_targets / weight_sum
  
  # Construction de la matrice bloc-diagonale optimale
  scaled_blocks_opt <- list()
  for (i in 1:n_blocks) {
    scaled_blocks_opt[[i]] <- block_alphas_opt[i] * xtx.raw.blocks[[i]]
  }
  xtx.raw.I_opt <- as.matrix(do.call(Matrix::bdiag, scaled_blocks_opt))
  
  # Construction de la matrice finale optimale
  xtx_opt <- xtx.raw.I_opt + alpha_non_diag_opt * xtx.raw.C + 
             (n_blocks - sum(block_alphas_opt)) * shrink.target_opt * diag(p)
  
  # Calcul du beta optimal
  beta.cov <- as.vector(scout::crossProdLasso(xtx_opt, xty, opt.lambda)$beta)
  beta.cov.lambda.max <- as.vector(scout::crossProdLasso(xtx_opt, xty, lmax)$beta)
  intercept <- ifelse(itcp, mean(y) - sum(beta.cov * (xm/xs)), 0)
  
  # Prédiction sur l'ensemble de test
  predict.test.values <- as.vector(x.test %*% beta.cov) + intercept
  
  # Calcul des métriques de performance
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
  
  # Préparation du résultat
  a <- list(
    'err' = DISCOM.tun.error, 
    'est.error' = DISCOM.est.error, 
    'lambda' = opt.lambda, 
    'alpha' = c(block_alphas_opt, alpha_non_diag_opt), 
    'train.error' = train.error, 
    'test.error' = DISCOM.test.error,
    'y.pred' = predict.test.values,
    'R2' = R2,
    'a0' = intercept, 
    'a1' = beta.cov, 
    'select' = select, 
    'xtx' = xtx_opt, 
    'fpr' = DISCOM.fpr, 
    'lambda.all' = lambda.all,
    'beta.cov.lambda.max' = beta.cov.lambda.max,
    'fnr' = DISCOM.fnr,
    "time" = as.numeric(time_taken[3])
  )
  
  return(a)
}
