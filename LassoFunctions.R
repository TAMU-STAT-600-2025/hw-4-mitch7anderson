# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  
  # [ToDo] Center and scale X
  Xmeans <- colMeans(X)
  Xcentered <- sweep(X, 2, Xmeans, FUN = "-")
  
  weights <- sqrt(colMeans(Xcentered^2))
  Xtilde <- sweep(Xcentered, 2, weights, FUN = "/")
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  s <- sign(a) * max(abs(a) - lambda, 0)
  
  return(s)
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  # Get sample size
  n <- nrow(Xtilde)
  r <- Ytilde - Xtilde %*% beta
  
  # Calculate objective
  fobj <- (1 / (2 * n)) * as.numeric(crossprod(r)) + lambda * sum(abs(beta))
 
  return(fobj)
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  if(nrow(Xtilde) != length(Ytilde)){
    stop("Xtilde and Ytilde must have same number of rows")
  }
  
  #[ToDo]  Check that lambda is non-negative
  if(lambda < 0){
    stop("Lambda must be non-negative")
  }
  
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  if(is.null(beta_start)){
    beta_start = rep(0, ncol(Xtilde))
  }else if(length(beta_start) != ncol(Xtilde)){
    stop("beta must be of length p")
  }
  
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  beta <- beta_start
  r <- Ytilde - drop(Xtilde %*% beta)             
  
  ## (1/n) * sum_j x_ij^2 for each column
  denom <- colSums(Xtilde * Xtilde) / n
  
  obj_old <- lasso(Xtilde, Ytilde, beta, lambda)
  
  repeat {
    for (j in seq_len(p)) {
      bj_old <- beta[j]
      # rho_j = (1/n) * x_j^T (r + x_j * b_j)  [temporarily add back j's contribution]
      rho_j <- as.numeric(crossprod(Xtilde[, j], r + Xtilde[, j] * bj_old)) / n
      
      # coordinate update
      bj_new <- soft(rho_j, lambda) / denom[j]
      
      if (bj_new != bj_old) {
        delta <- bj_new - bj_old
        r <- r - Xtilde[, j] * delta               
        beta[j] <- bj_new
      }
    }
    
    fmin <- lasso(Xtilde, Ytilde, beta, lambda)
    # Check stoping condition
    if ((obj_old - fmin) < eps) break           
    obj_old <- fmin
  }
  
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  if(nrow(Xtilde) != length(Ytilde)){
    stop("Xtilde and Ytilde must have same number of rows")
  }
 
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  if (!is.null(lambda_seq)) {
    lambda_seq <- sort(unique(lambda_seq[is.finite(lambda_seq) & lambda_seq >= 0]), decreasing = TRUE)
    if (length(lambda_seq) == 0L) {
      warning("No valid nonnegative values in supplied lambda_seq; proceeding as if values were not supplied.")
      lambda_seq <- NULL
    }
  }
  
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
  if (is.null(lambda_seq)) {
    n <- nrow(Xtilde)
    lambda_max <- max(abs(as.numeric(crossprod(Xtilde, Ytilde)) / n))
    if (is.finite(lambda_max) && lambda_max > 0) {
      lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda))
    } else {
      # fallback if lambda_max == 0 or not finite
      lambda_seq <- rep(0, n_lambda)
    }
    # ensure largest -> smallest
    lambda_seq <- sort(lambda_seq, decreasing = TRUE)
  }
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  p <- ncol(Xtilde)
  L <- length(lambda_seq)
  beta_mat <- matrix(0, nrow = p, ncol = L)
  fmin_vec <- numeric(L)
  
  beta_start <- rep(0, p)
  for (i in seq_len(L)) {
    lam_i <- lambda_seq[i]
    fit_i <- fitLASSOstandardized(Xtilde, Ytilde, lambda = lam_i, beta_start = beta_start, eps = eps)
    beta_mat[ , i] <- fit_i$beta
    fmin_vec[i]   <- fit_i$fmin
    # Warm start for next lambda
    beta_start    <- fit_i$beta 
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  std <- standardizeXY(X, Y)
  Xtilde  <- std$Xtilde
  Ytilde  <- std$Ytilde
  Ymean   <- std$Ymean
  Xmeans  <- std$Xmeans
  weights <- std$weights
 
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  path <- fitLASSOstandardized_seq(Xtilde, Ytilde,
                                   lambda_seq = lambda_seq,
                                   n_lambda   = n_lambda,
                                   eps        = eps)
 
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  lambda_seq <- path$lambda_seq
  beta_std_mat <- path$beta_mat   # p x L in standardized scale
  L <- length(lambda_seq)
  p <- ncol(X)
  
  beta_mat <- matrix(0, nrow = p, ncol = L)
  nz_w <- weights != 0
  if (any(nz_w)) {
    beta_mat[nz_w, ] <- sweep(beta_std_mat[nz_w, , drop = FALSE], 1, weights[nz_w], "/")
  }
  
  beta0_vec <- as.numeric(Ymean - crossprod(Xmeans, beta_mat))
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  # [ToDo] Fit Lasso on original data using fitLASSO
  n <- nrow(X)
  p <- ncol(X)
  
  full_fit <- fitLASSO(X, Y, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  lambda_seq <- full_fit$lambda_seq      
  L <- length(lambda_seq)
  beta_mat <- full_fit$beta_mat          
  beta0_vec <- full_fit$beta0_vec
 
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  if (is.null(fold_ids)) {
    # balanced random folds 1..k
    fold_ids <- sample(rep(1:k, length.out = n))
  } else {
    # coerce to integer
    fold_ids <- as.integer(fold_ids)
    uniq <- sort(unique(fold_ids))
    # map arbitrary labels to 1..K
    fold_ids <- match(fold_ids, uniq)
    k <- length(uniq)
  }
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  mse_mat <- matrix(NA_real_, nrow = k, ncol = L)
  
  for (fold in seq_len(k)) {
    idx_val <- which(fold_ids == fold)
    idx_tr  <- setdiff(seq_len(n), idx_val)
    
    Xtr <- X[idx_tr, , drop = FALSE]
    Ytr <- Y[idx_tr]
    Xvl <- X[idx_val, , drop = FALSE]
    Yvl <- Y[idx_val]
    
    # Fit train split with THE SAME lambda_seq (and same eps)
    fit_tr <- fitLASSO(Xtr, Ytr, lambda_seq = lambda_seq, n_lambda = L, eps = eps)
    
    # Predict on validation split for every lambda: yhat = beta0 + Xvl %*% beta
    # (fitLASSO already returns original-scale betas/intercepts, so just use them)
    preds <- sweep(Xvl %*% fit_tr$beta_mat, 2, fit_tr$beta0_vec, FUN = "+")  # n_val x L
    # MSE per lambda
    resid <- preds - Yvl
    mse_fold <- colMeans(resid^2)
    mse_mat[fold, ] <- mse_fold
  }
  
  # CV summary
  cvm  <- colMeans(mse_mat)                            # mean MSE across folds
  cvse <- apply(mse_mat, 2, sd) / sqrt(k)              # standard error across folds
  
  # [ToDo] Find lambda_min
  i_min <- which.min(cvm)
  lambda_min <- lambda_seq[i_min]

  # [ToDo] Find lambda_1SE
  thresh <- cvm[i_min] + cvse[i_min]
  # choose the largest lambda with cvm <= thresh
  eligible <- which(cvm <= thresh)
  if (length(eligible)) {
    i_1se <- max(eligible)
  } else {
    i_1se <- i_min
  }
  lambda_1se <- lambda_seq[i_1se]
  
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

