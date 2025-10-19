source("LassoFunctions.R")

n = 50 # sample size
beta = c(1, 0.5) # true coefficients vector
beta0 = 2 # intercept, nonzero
p = length(beta)
sigma = 0.4 # noise standard deviation

library(mnormt)
set.seed(983645) # set seed

# Generate matrix of covariates
Sigma = matrix(0.7, p, p) + diag(rep(1-0.7, p)) # covariance for design X
X = rmnorm(n, mean = rep(0, p), varcov = Sigma) # matrix of covariates

# Generate response Y
Y = beta0 + X %*% beta + sigma * rnorm(n) 

std = standardizeXY(X, Y)
crossprod(std$Xtilde[ , 1])/n

## ---------------------------
## Helpers
## ---------------------------
std_XY <- function(X, Y) {
  Xc <- scale(X, center = TRUE, scale = TRUE)
  Yc <- as.numeric(scale(Y, center = TRUE, scale = FALSE))
  list(X = Xc, Y = Yc)
}

# KKT residual for LASSO (standardized, no intercept)
# If beta_j == 0: |(1/n) x_j^T r| <= lambda
# If beta_j != 0: (1/n) x_j^T r == lambda * sign(beta_j)
kkt_violation <- function(X, Y, beta, lambda, tol = 1e-6) {
  n <- nrow(X)
  r <- Y - drop(X %*% beta)
  g <- as.numeric(crossprod(X, r)) / n # (1/n) X^T r
  nz <- which(beta != 0)
  z  <- which(beta == 0)
  v1 <- if (length(nz)) max(abs(g[nz] - lambda * sign(beta[nz]))) else 0
  v2 <- if (length(z))  max(pmax(0, abs(g[z]) - lambda)) else 0
  max(v1, v2, na.rm = TRUE)
}

## ---------------------------
## T1: Large lambda => beta = 0 (KKT at zero)
## ---------------------------
set.seed(1)
n <- 80; p <- 30
X <- matrix(rnorm(n * p), n, p)
beta_true <- rep(0, p); beta_true[1:5] <- rnorm(5, 0, 2)
Y <- X %*% beta_true + rnorm(n, 0, 1.0)

s <- std_XY(X, Y); Xs <- s$X; Ys <- s$Y
n <- nrow(Xs)

# KKT at beta=0 says choose lambda >= max(|(1/n) x_j^T y|)
lam_big <- max(abs(as.numeric(crossprod(Xs, Ys)) / n)) + 1e-4

fit1 <- fitLASSOstandardized(Xs, Ys, lambda = lam_big, eps = 1e-8)
stopifnot(all.equal(unname(fit1$beta), rep(0, p), tolerance = 1e-10))

# KKT residual should be ~0
viol1 <- kkt_violation(Xs, Ys, fit1$beta, lam_big)
stopifnot(viol1 < 1e-6)

## ---------------------------
## T2: lambda = 0 => OLS (no intercept; X, Y centered)
## ---------------------------
lam0 <- 0
fit2 <- fitLASSOstandardized(Xs, Ys, lambda = lam0, eps = 1e-10)

# Compare predictions to least-squares via QR
qrX <- qr(Xs)
beta_ls <- qr.coef(qrX, Ys)       # least-squares solution
beta_ls[is.na(beta_ls)] <- 0      # handle p>n minimal-norm via qr
pred_fit2 <- drop(Xs %*% fit2$beta)
pred_ls   <- drop(Xs %*% beta_ls)
stopifnot(all.equal(pred_fit2, pred_ls, tolerance = 1e-6))

## ---------------------------
## T3: Orthogonal design => closed-form soft-threshold solution
## Build X with orthonormal columns and scale to ||x_j||^2 = n (=> denom_j = 1)
## ---------------------------
set.seed(2)
n <- 60; p <- 20
Z <- matrix(rnorm(n * p), n, p)
Q <- qr.Q(qr(Z))                 # n x n orthonormal
Xo <- Q[, 1:p, drop = FALSE]     # n x p with orthonormal columns
Xo <- Xo * sqrt(n)               # now each col has sum(x^2) = n -> denom = 1
Yo <- rnorm(n)

lambda <- 0.25
# Closed-form: beta_j = soft((1/n) x_j^T y, lambda) / denom_j, but denom_j = 1 here
z <- as.numeric(crossprod(Xo, Yo)) / n
beta_cf <- soft(z, lambda)

# Our solver
fit3 <- fitLASSOstandardized(Xo, Yo, lambda = lambda, eps = 1e-10)
stopifnot(all.equal(unname(fit3$beta), beta_cf, tolerance = 1e-7))

## ---------------------------
## T4: Convergence / monotonicity
## Re-run with warm start; objective should not drop by >= eps beyond the first solution
## ---------------------------
obj1 <- lasso(Xs, Ys, fit1$beta, lam_big)
fit1b <- fitLASSOstandardized(Xs, Ys, lambda = lam_big, beta_start = fit1$beta, eps = 1e-8)
obj1b <- lasso(Xs, Ys, fit1b$beta, lam_big)
stopifnot((obj1 - obj1b) < 1e-8 + 1e-12)

## ---------------------------
## T5: p >> n stress (sparsity + KKT)
## ---------------------------
set.seed(3)
n <- 100; p <- 5000
Xbig <- matrix(rnorm(n * p), n, p)
beta_true <- rep(0, p); beta_true[1:10] <- rnorm(10, 0, 2)
Ybig <- Xbig %*% beta_true + rnorm(n, 0, 1)

s <- std_XY(Xbig, Ybig); Xbs <- s$X; Ybs <- s$Y
n <- nrow(Xbs)

# Choose lambda to yield a sparse solution (heuristic)
lam <- quantile(abs(as.numeric(crossprod(Xbs, Ybs)) / n), 0.9)

fit5 <- fitLASSOstandardized(Xbs, Ybs, lambda = lam, eps = 1e-4)
# KKT residual should be small
viol5 <- kkt_violation(Xbs, Ybs, fit5$beta, lam, tol = 1e-6)
stopifnot(viol5 < 1e-3)
# Sanity: solution is sparse
stopifnot(sum(fit5$beta != 0) < n)   # usually << n

cat("All tests passed ✅\n")

## ---------------------------
## T6 (optional): glmnet cross-check if available
## ---------------------------
if (requireNamespace("glmnet", quietly = TRUE)) {
  library(glmnet)
  set.seed(4)
  n <- 120; p <- 60
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- rep(0, p); beta_true[1:8] <- rnorm(8, 0, 1.5)
  Y <- X %*% beta_true + rnorm(n, 0, 1)
  
  s <- std_XY(X, Y); Xs <- s$X; Ys <- s$Y
  n <- nrow(Xs)
  
  # Match our objective scaling: glmnet uses (1/(2n))||r||^2 by default if we set standardize=FALSE and no intercept
  lam <- 0.2
  fit_me <- fitLASSOstandardized(Xs, Ys, lambda = lam, eps = 1e-8)
  
  fit_glm <- glmnet::glmnet(Xs, Ys, alpha = 1, lambda = lam, intercept = FALSE, standardize = FALSE)
  beta_glm <- as.numeric(fit_glm$beta)
  # Compare predictions (coefs can differ when features are highly correlated, but preds should be close)
  pred_me  <- drop(Xs %*% fit_me$beta)
  pred_glm <- drop(Xs %*% beta_glm)
  stopifnot(all.equal(pred_me, pred_glm, tolerance = 1e-5))
  
  cat("glmnet cross-check passed ✅\n")
}

gen_problem <- function(n = 100, p = 30, s = 6, snr = 4, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- numeric(p); beta_true[seq_len(s)] <- rnorm(s, 0, 1.5)
  f <- drop(X %*% beta_true)
  sigma <- sqrt(var(f) / snr)
  Y <- f + rnorm(n, 0, sigma)
  list(X = X, Y = Y, beta_true = beta_true)
}

prob <- gen_problem(n = 100, p = 30, s = 6, seed = 11)
prob$X
std = standardizeXY(prob$X, prob$Y)

colSums(std$Xtilde * std$Xtilde) / n

## ==========================================
## T1: Default path — shapes, order, λ_max zero-soln, KKT mid-path
## ==========================================
t1 <- function() {
  cat("T1: default path sanity… ")
  prob <- gen_problem(n = 100, p = 30, s = 6, seed = 11)
  std  <- standardizeXY(prob$X, prob$Y)
  
  out <- fitLASSOstandardized_seq(std$Xtilde, std$Ytilde,
                                  lambda_seq = NULL, n_lambda = 50, eps = 1e-8)
  
  # shapes
  stopifnot(length(out$lambda_seq) == 50)
  stopifnot(identical(dim(out$beta_mat), c(ncol(prob$X), 50)))
  stopifnot(length(out$fmin_vec) == 50)
  
  # λ descending
  stopifnot(isTRUE(all(diff(out$lambda_seq) <= 0)))
  
  # first β should be zero at λ_max (up to numerical noise)
  stopifnot(max(abs(out$beta_mat[, 1])) < 1e-10)
  
  # KKT small at a mid λ
  k <- 20
  viol <- kkt_violation(std$Xtilde, std$Ytilde, out$beta_mat[, k], out$lambda_seq[k])
  stopifnot(viol < 1e-5)
  
  cat("OK\n")
}

t1()


