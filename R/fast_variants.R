#' Fast AdapDiscom
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
#' @param robust Integer, 0 for classical, 1 for robust estimation
#' @param n.l Integer, number of tuning parameter (`l`) values for fast variants
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
fast_adapdiscom <- function(beta, x, y, x.tuning, y.tuning, x.test, y.test, nlambda, nalpha, pp,
                            robust = 0, n.l = 30, standardize = TRUE, itcp = TRUE,
                            lambda.min.ratio = NULL,  k.value = 1.5) {
  
  start_time <- proc.time()
  
  x_std = standardize_x(as.matrix(x), robust = robust, k.value = k.value)
  xm <- x_std$x.mean
  xs <- x_std$x.sd
  if (standardize) {
    x <- x_std$x
    x.tuning <- fit_standardize_x(x.tuning, xm, xs)
    x.test <- fit_standardize_x(x.test, xm, xs)
  }
  
  
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  n.tuning <- dim(x.tuning)[1]
  n.test <- dim(x.test)[1]
  n_blocks <- length(pp)
  
  alpha.all <- 10^seq(10^(-10), -3, length = nalpha)
  
  lambda.max <- lambda_max(x, y, Methode = "discom", robust = robust)
  nobs <- dim(na.omit(x))[1]
  if(is.null(lambda.min.ratio) || length(lambda.min.ratio) == 0 || lambda.min.ratio <= 0 ||
     lambda.min.ratio > 1){
    lambda.min.ratio <- ifelse(nobs < p, 0.01, 1e-04)
  }
  
  lambda.min <- lambda.max * lambda.min.ratio
  lambda.all <- exp(seq(log(lambda.max), log(lambda.min), length.out = nlambda))
  
  xtx.raw <- compute.xtx(x, robust = robust, k_value = k.value)
  xty <- compute.xty(x, y, robust = robust, k_value = k.value)
  
  # Création des indices pour les sous-matrices
  indices <- get_block_indices(pp)
  
  # Création des sous-matrices individuelles
  xtx.raw.blocks <- list()
  shrink.targets <- numeric(n_blocks)
  nj <- numeric(n_blocks)
  x_dummy = x
  # Replace Non NA by 1
  x_dummy[!is.na(x_dummy)] <- 1
  # Replace NA by 0
  x_dummy[is.na(x_dummy)] <- 0
  
  njt <- min(x_dummy%*%t(x_dummy))
  
  for (i in 1:n_blocks) {
    idx_range <- indices$starts[i]:indices$ends[i]
    xtx.raw.blocks[[i]] <- xtx.raw[idx_range, idx_range]
    shrink.targets[i] <- sum(diag(xtx.raw.blocks[[i]]))/p
    nj[i] <-  min(colSums(!is.na(x[, idx_range])))
  }
  
  mk <- log(p)/nj
  mk <- sqrt(mk)
  mc <- sqrt(log(p)/njt)
  
  # Création de la matrice bloc-diagonale initiale
  xtx.raw.I <- as.matrix(do.call(Matrix::bdiag, xtx.raw.blocks))
  
  # Calcul de xtx.raw.C
  xtx.raw.C <- xtx.raw - xtx.raw.I
  
  # Calcul du shrink.target global
  shrink.target <- sum(diag(xtx.raw))/p
  
  # Calcul des paramètres de fast adapdiscom
  lmax <- 1/mc
  Mat_1 <- Map(function(x, y) {
    x * y
  }, x = xtx.raw.blocks, y = as.list(mc - mk))
  Mat_1 <- do.call(Matrix::bdiag, Mat_1)
  
  sum_square_mk <- sum(mk^2)
  Mat_2 <- Map(function(x, y) {
    (sum(diag(x)) * (y^2)) / sum_square_mk
  }, x = xtx.raw.blocks, y = as.list(mk))
  Mat_2 <- sum(unlist(Mat_2))
  Mat_2 <- sum(mk * Mat_2)
  Mat_2 <- (Mat_2 * diag(p))/p
  
  Mat <- Mat_1 + Mat_2
  
  kappa_sigma <- min(eigen(xtx.raw, only.values = TRUE)$values)
  kappa_mat <- min(eigen(Mat, only.values = TRUE)$values)
  
  lmin <- (kappa_sigma)/(mc*kappa_sigma - kappa_mat)
  lmin <- max(lmin, 0)
  l.all <- seq(lmin, lmax, length = n.l)
  
  dim_array <- c(n.l, nlambda)
  DISCOM.tun.error <- array(NA, dim = dim_array)
  
  for (row in 1:n.l) {
    l0 <- l.all[row]
    current_alphas <- 1 - l0 * mk
    
    # Les n_blocks premières alphas sont pour les blocs diagonaux
    block_alphas <- current_alphas
    
    # Le dernier alpha est pour xtx.raw.C
    alpha_non_diag <- 1 - l0 * mc
    
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
    
    # Supprimer scaled_blocks pour libérer de la mémoire
    rm(scaled_blocks)
    
    # Construction de la matrice finale
    xtx <- xtx.raw.I_scaled + alpha_non_diag * xtx.raw.C + 
           (n_blocks - sum(block_alphas)) * shrink.target_weighted * diag(p)
    
    # Vérification de la positivité
    if (any(is.na(xtx)) || any(is.infinite(as.vector(xtx)))) {
      DISCOM.tun.error[row, ] <- 10^8
    } else {
      eigen_values <- try(eigen(xtx, only.values = TRUE)$values, silent = TRUE)
      if (inherits(eigen_values, "try-error") || min(eigen_values) < 0) {
        DISCOM.tun.error[row, ] <- 10^8
      } else {
        beta.initial <- rep(0, p)
        for (k in 1:nlambda) {
          beta.cov <- as.vector(scout::crossProdLasso(xtx, xty, lambda.all[k], 
                                                       beta.init = beta.initial)$beta)
          beta.initial <- beta.cov
          intercept <- ifelse(itcp, mean(y) - sum(beta.cov * (xm/xs)), 0)
          DISCOM.tun.values <- as.vector(as.matrix(x.tuning) %*% beta.cov) + intercept
          mse <- mean((y.tuning - DISCOM.tun.values)^2)
          DISCOM.tun.error[row, k] <- mse
        }
      }
    }
  }
  
  # Trouver les paramètres optimaux
  opt.index <- as.vector(which(DISCOM.tun.error == min(DISCOM.tun.error), arr.ind = TRUE)[1, ])
  train.error <- min(DISCOM.tun.error)
  
  # Extraire les valeurs alpha optimales
  opt.alphas <- 1 - l.all[opt.index[1]] * mk
  opt.lambda <- lambda.all[opt.index[2]]
  
  # Extraire les valeurs alpha optimales pour les blocs
  block_alphas_opt <- opt.alphas[1:n_blocks]
  alpha_non_diag_opt <- 1 - l.all[opt.index[1]] * mc
  
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
  beta.cov.lambda.max <- as.vector(scout::crossProdLasso(xtx_opt, xty, lambda.max)$beta)
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

#' Fast DISCOM
#'
#' @param beta Vector, true beta coefficients (optional)
#' @param x Matrix, training data
#' @param y Vector, training response
#' @param x.tuning Matrix, tuning data
#' @param y.tuning Vector, tuning response
#' @param x.test Matrix, test data
#' @param y.test Vector, test response
#' @param nlambda Integer, number of lambda values
#' @param pp Vector, block sizes
#' @param robust Integer, 0 for classical, 1 for robust estimation
#' @param n.l Integer, number of tuning parameter (`l`) values for fast variants
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
fast_discom <- function(beta, x, y, x.tuning, y.tuning, x.test, y.test, nlambda, pp,
                        robust = 0, n.l = 30, standardize = TRUE, itcp = TRUE,
                        lambda.min.ratio = NULL,  k.value = 1.5) {
  
  start_time <- proc.time()
  x_std = standardize_x(as.matrix(x), robust = robust, k.value = k.value)
  xm <- x_std$x.mean
  xs <- x_std$x.sd
  if (standardize) {
    x <- x_std$x
    x.tuning <- fit_standardize_x(x.tuning, xm, xs)
    x.test <- fit_standardize_x(x.test, xm, xs)
  }
  
  # Récupérer les dimensions
  p <- ncol(x)
  p1 <- pp[1]
  p2 <- pp[2]
  if (length(pp) > 2) p3 <- pp[3]
  if (length(pp) > 3) p4 <- pp[4]
  
  beta.true <- beta
  
  # Calcul des matrices croisées
  xtx.raw <- compute.xtx(x, robust = robust, k_value = k.value)
  xty <- compute.xty(x, y, robust = robust, k_value = k.value)
  
  if (length(pp) == 4) {
    xtx.raw.I <- as.matrix(Matrix::bdiag(xtx.raw[1:p1, 1:p1], 
                                        xtx.raw[(p1+1):(p1+p2), (p1+1):(p1+p2)],
                                        xtx.raw[(p1+p2+1):(p1+p2+p3), (p1+p2+1):(p1+p2+p3)],
                                        xtx.raw[(p1+p2+p3+1):(p1+p2+p3+p4), (p1+p2+p3+1):(p1+p2+p3+p4)]))
    xtx.raw.C <- xtx.raw - xtx.raw.I
    shrink.target <- sum(diag(xtx.raw))/p
  } else if (length(pp) == 3) {
    xtx.raw.I <- as.matrix(Matrix::bdiag(xtx.raw[1:p1, 1:p1], 
                                        xtx.raw[(p1+1):(p1+p2), (p1+1):(p1+p2)],
                                        xtx.raw[(p1+p2+1):(p1+p2+p3), (p1+p2+1):(p1+p2+p3)]))
    xtx.raw.C <- xtx.raw - xtx.raw.I
    shrink.target <- sum(diag(xtx.raw))/p
  } else if (length(pp) == 2) {
    xtx.raw.I <- as.matrix(Matrix::bdiag(xtx.raw[1:p1, 1:p1], 
                                        xtx.raw[(p1+1):(p1+p2), (p1+1):(p1+p2)]))
    xtx.raw.C <- xtx.raw - xtx.raw.I
    shrink.target <- sum(diag(xtx.raw))/p
  }
  
  nj <- min(colSums(!is.na(x)))
  x_dummy = x
  # Replace Non NA by 1
  x_dummy[!is.na(x_dummy)] <- 1
  # Replace NA by 0
  x_dummy[is.na(x_dummy)] <- 0
  
  njt <- min(x_dummy%*%t(x_dummy))
  rm(x_dummy)
  # Définir les constantes
  c1 <- sqrt(log(p)/nj)
  c2 <- sqrt(log(p)/njt)
  
  # Créer la matrice A
  A.matrix <- (c2-c1)*xtx.raw.I + c1*shrink.target*diag(p)
  
  # Calculer les limites pour K
  K.min <- max(0, min(eigen(xtx.raw)$values)/(c2*min(eigen(xtx.raw)$values) - min(eigen(A.matrix)$values)))
  K.max <- 1/c2
  
  # Générer les séquences pour K et lambda
  K.all <- seq(K.min, K.max, length = n.l)
  n <- dim(x)[1]
  lambda.max <- lambda_max(x, y, Methode = "discom", robust = robust)
  
  nobs <- dim(na.omit(x))[1]
  if(is.null(lambda.min.ratio) || length(lambda.min.ratio) == 0 || lambda.min.ratio <= 0 ||
     lambda.min.ratio > 1){
    lambda.min.ratio <- ifelse(nobs < p, 0.01, 1e-04)
  }
  
  lambda.min <- lambda.max * lambda.min.ratio
  lambda.all <- exp(seq(log(lambda.max), log(lambda.min), length.out = nlambda))
  
  # Matrice pour stocker les erreurs de tuning
  DISCOM.tun.error <- matrix(NA, n.l, nlambda)
  
  # Boucle principale pour différentes valeurs de K et lambda
  for (i in 1:n.l) {
    alpha1 <- 1 - K.all[i]*c1
    alpha2 <- 1 - K.all[i]*c2
    xtx <- alpha1*xtx.raw.I + alpha2*xtx.raw.C + (1-alpha1)*shrink.target*diag(p)
    
    beta.initial <- rep(0, p)
    
    for (k in 1:nlambda) {
      beta.cov <- as.vector(scout::crossProdLasso(xtx, xty, lambda.all[k], 
                                                       beta.init = beta.initial)$beta)
      beta.initial <- beta.cov
      intercept <- ifelse(itcp, mean(y) - sum(beta.cov * (xm/xs)), 0)
      
      # Calculer les valeurs prédites et l'erreur
      DISCOM.tun.values <- as.vector(as.matrix(x.tuning) %*% beta.cov) + intercept
      DISCOM.tun.error[i, k] <- mean((y.tuning - DISCOM.tun.values)^2)
    }
  }
  
  # Trouver l'index optimal
  opt.index <- as.vector(which(DISCOM.tun.error == min(DISCOM.tun.error), arr.ind = TRUE)[1, ])
  
  # Calculer les alphas optimaux
  opt.alpha1 <- 1 - K.all[opt.index[1]]*c1
  opt.alpha2 <- 1 - K.all[opt.index[1]]*c2
  opt.lambda <- lambda.all[opt.index[2]]
  
  # Calculer la matrice xtx optimale
  xtx_opt <- opt.alpha1*xtx.raw.I + opt.alpha2*xtx.raw.C + (1-opt.alpha1)*shrink.target*diag(p)
  
  # Calculer le beta final
  beta.cov <- as.vector(scout::crossProdLasso(xtx_opt, xty, opt.lambda)$beta)
  beta.cov.lambda.max <- as.vector(scout::crossProdLasso(xtx_opt, xty, lambda.max)$beta)
  intercept <- ifelse(itcp, mean(y) - sum(beta.cov * (xm/xs)), 0)
  
  # Prédire sur l'ensemble de test
  predict.test.values <- as.vector(x.test %*% beta.cov) + intercept
  
  # Calculer l'erreur de test
  DISCOM.test.error <- mean((y.test - predict.test.values)^2)
  
  # Calculer les métriques de performance
  if (!is.null(beta)) {
    DISCOM.fpr <- sum((beta == 0) & (beta.cov != 0))/sum(beta == 0)
    DISCOM.fnr <- sum((beta != 0) & (beta.cov == 0))/sum(beta != 0)
    DISCOM.est.error <- sqrt(sum((beta.cov - beta)^2))
  } else {
    DISCOM.fpr <- NA
    DISCOM.fnr <- NA
    DISCOM.est.error <- NA
  }
  
  # Calculer l'erreur d'entraînement
  train.error <- min(DISCOM.tun.error)
  
  # Calculer R²
  R2 <- cor(predict.test.values, y.test)^2
  
  # Déterminer les variables sélectionnées
  select <- sum(as.vector(as.integer(beta.cov != 0)))
  
  # Calculer les alphas pour les blocs
  block_alphas_opt <- c(opt.alpha1, opt.alpha2)
  alpha_non_diag_opt <- opt.alpha2
  
  end_time <- proc.time()
  time_taken <- end_time - start_time
  
  # Créer la liste de résultats
  results <- list(
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
    'lambda.all' = lambda.all,
    'beta.cov.lambda.max' = beta.cov.lambda.max,
    'xtx' = xtx_opt, 
    'fpr' = DISCOM.fpr, 
    'fnr' = DISCOM.fnr,
    "time" = as.numeric(time_taken[3])
  )
  
  return(results)
}