################################################################################
################################## dsplineML ###################################
################################################################################
#' @importFrom nloptr nloptr
#' @importFrom numDeriv jacobian
#' @importFrom splines splineDesign

dsplineML <- function(X, n, k, boundary_effects = TRUE, only_theta = FALSE, init_int.knt = NULL,
                      init_theta = NULL, optm = "local") {

  # Initial checks
  if (only_theta && is.null(init_int.knt)) stop("For ML estimation of coeffients a knot vector should be provided.")
  if (!is.null(init_int.knt) && length(init_int.knt) != k){
    k <- length(init_int.knt)
    warning("k corresponds to the # of internal knots.")
  }
  if (k == 0) only_theta <- TRUE

  p <- n + k
  extr <- range(X)

  # Define f_X, the density function
  f_X <- function(params) {
    p <- n + k
    theta <- params[1:p]
    if(!only_theta) int.knt <- params[(p+1):(p+k)]
    knt <- sort(c(rep(extr,n), int.knt))

    # Compute the basis matrix using current knots
    basisMatrix <- splineDesign(knots = knt, x = X, ord = n, derivs = rep(0, length(X)), outer.ok = TRUE)

    return(basisMatrix %*% theta)  # Compute f(X) = basisMatrix * theta
  }

  objFun <- function(params) {
    theta <- params[1:p]
    if(!only_theta) int.knt <- params[(p + 1):(p + k)]

    # Compute basis matrix and f_values
    knt <- sort(c(rep(extr, n), int.knt))
    basisMatrix <- splineDesign(knots = knt, x = X, ord = n, derivs = rep(0, length(X)), outer.ok = TRUE)
    f_values <- pmax(basisMatrix %*% theta, 1e-3)

    # Log-likelihood computation
    logLik <- sum(log(f_values))

    # Gradient w.r.t. theta
    grad_theta <- - colSums( basisMatrix / as.vector(f_values) )
    # Gradient w.r.t. knots
    # grad_knt <- rep(0, k)
    # for (i in 1:k) {
    #   # Compute derivative of basis functions w.r.t. the i-th knot
    #   deriv_basis <- splineDesign(knots = knt, x = X, ord = n, derivs = rep(1, length(X)), outer.ok = TRUE)
    #   grad_knt[i] <- sum((1 / f_values) * (deriv_basis %*% theta))
    # }

    if(only_theta) {
      grad_knt <- NULL
    } else {
      grad_knt <- - t(jacobian(func = function(knt) f_X(c(theta, knt)), x = int.knt)) %*% (1 / f_values)
    }

    # Combine gradients
    grad <- c(grad_theta, grad_knt)

    return(list("objective" = -logLik, "gradient" = grad))
  }

  # 2. Define constraints
  # 2.1. Inequality constraints
  ineq_constraints <- function(params) {
    # re-formulate ineq constraints to be of form g(x) <= 0
    theta <- params[1:p]
    if(!only_theta) int.knt <- params[(p + 1):(p + k)]

    # (1) θ_j ≥ 0, j = 1,...,p
    constr1 <- -theta
    # Initialize the matrix with zeros
    grad1 <- matrix(0, nrow = p, ncol = length(params))
    # Use matrix indexing to fill in 1s and -1s
    for (i in 1:p) {
      grad1[i, i] <- -1
    }

    if(only_theta){
      constr2 <- constr3 <- grad2 <- grad3 <- NULL
    } else {

      # (2) t_{n+1}<...<t_{n+k}
      if(length(int.knt) > 1) {
        epsilon <- 1e-6 # to ensure strict inequality
        constr2 <- c(epsilon - diff(int.knt))
        # Initialize the matrix with zeros
        grad2 <- matrix(0, nrow = k-1, ncol = p+k)
        for (i in 1:(k-1)) {
          grad2[i, p + i] <- 1
          grad2[i, p + i + 1] <- -1
        }
      } else {
        constr2 <- grad2 <- NULL
      }

      # (3) t_{1}=...=t_{n} < t_{n+1} & t_{n+k} < t_{n+k+1}=...=t_{2n+k}
      epsilon <- 1e-6 # to ensure strict inequality
      diff1 <- int.knt[1]-extr[1]
      diff2 <- extr[2]-int.knt[k]
      constr3 <- c(epsilon - diff1,
                   epsilon - diff2)
      # Initialize the matrix with zeros
      grad3 <- matrix(0, nrow = 2, ncol = p+k)
      grad3[1, p+1] <- -1
      grad3[2, p+k] <- 1

    }

    return( list( "constraints" = c(constr1, constr2, constr3), "jacobian" = rbind(grad1, grad2, grad3) ) )
  }

  # 2.2. Equality constraints
  eq_constraints <- function(params) {
    theta <- params[1:p]
    if(!only_theta) int.knt <- params[(p + 1):(p + k)]
    knt <- sort(c(rep(extr, n), int.knt))

    # (1) ∑ θ_j * (t_{j+n} - t_j) / n = 1
    constr1 <- sum(theta * diff(knt, n)/n) - 1

    # Gradient w.r.t. theta
    grad_theta <- diff(knt, n)/n
    # Gradient w.r.t. knots
    if (only_theta) {
      grad_knt <- NULL
    } else {
      grad_knt <- (-theta/n)[(n+1):(n+k)]
    }

    # Initialize the matrix with zeros
    grad1 <- c(grad_theta, grad_knt)

    if (boundary_effects) {
      # (2) θ_1 = θ_p = 0
      constr2 <- c(theta[1] - 0,
                   theta[p] - 0)
      # Initialize the matrix with zeros
      grad2 <- matrix(0, nrow = 2, ncol = length(params))
      # Set the first element of the first row to 1
      grad2[1, 1] <- 1
      # Set the last element of the last row to 1
      grad2[2, p] <- 1
    } else {
      constr2 <- grad2 <- NULL
    }

    return( list( "constraints" = c(constr1, constr2), "jacobian" = rbind(grad1, grad2) ) )
  }



  # Initialize knots (quantile-based placement)
  if (only_theta) {
    int.knt <- init_int.knt; init_int.knt <- NULL
    init_knt <- sort(c(rep(extr,n), int.knt))
  } else {
    if(is.null(init_int.knt)) {
      init_int.knt <- quantile(X, probs = seq(0, 1, length.out = k + 2))[-c(1, k + 2)]
      init_int.knt <- sort(unique(init_int.knt))  # Ensure sorted and unique
      init_knt <- sort(c(rep(extr,n), init_int.knt))
    } else {
      init_knt <- sort(c(rep(extr,n), init_int.knt))
    }
  }
  # Initialize coefficients
  if(is.null(init_theta)) {
    if (boundary_effects) {
      init_theta <-  c(1e-3, rep(n / sum(diff(init_knt, n)[-c(1, p)]), p - 2), 1e-3)
    } else {
      init_theta <- rep(n / sum(diff(init_knt, n)), p)
    }
  }

  # Combine into parameter vector
  init_params <- c(init_theta, init_int.knt)


  # # Perform Constrained Optimization using nloptr
  if (optm == "global") {

    # Global optimizer
    opts_global <- list(
      "algorithm" = "NLOPT_GN_ISRES",
      "xtol_rel" = 1.0e-12,
      "maxeval" = 10000
    )
    opt_result <- nloptr(
      x0 = init_params,
      eval_f = objFun,
      eval_g_ineq = ineq_constraints,
      eval_g_eq = eq_constraints,
      opts = opts_global
    )
  } else if (optm == "local") {

    local_opts <- list(
      "algorithm" = "NLOPT_LD_MMA",
      "xtol_rel" = 1.0e-10
    )
    opts <- list(
      "algorithm" = "NLOPT_LD_AUGLAG",
      "xtol_rel" = 1.0e-10,
      "maxeval" = 10000,
      "local_opts" = local_opts
    )

    # # Local optimizer
    # opts_local <- list(
    #   "algorithm" = "NLOPT_LD_SLSQP",
    #   "xtol_rel" = 1.0e-12,
    #   "maxeval" = 10000
    # )

    # local_opts <- list(
    #   "algorithm" = "NLOPT_LD_LBFGS",
    #   "xtol_rel" = 1.0e-10
    # )
    #
    # opts <- list(
    #   "algorithm" = "NLOPT_LD_AUGLAG",
    #   "xtol_rel" = 1.0e-10,
    #   "maxeval" = 10000,
    #   "local_opts" = local_opts
    # )

    opt_result <- nloptr(
      x0 = init_params,
      eval_f = objFun,
      eval_g_ineq = ineq_constraints,
      eval_g_eq = eq_constraints,
      opts = opts
    )

  }

  # Extract results
  theta_hat <- opt_result$solution[1:p]
  if(!only_theta) int.knt <- opt_result$solution[(p+1):(p+k)]
  knt_hat <- sort(c(rep(extr,n),int.knt))

  basisMatrix <- splineDesign(knots = knt_hat, x = X, ord = n, derivs = rep(0, length(X)), outer.ok = TRUE)
  f_X_hat <- basisMatrix %*% theta_hat

  # Check
  if( all( all(theta_hat >= 0) &&
           abs(sum(theta_hat * diff(knt_hat, n)/n) - 1) < 1e-9  &&
           all( diff(knt_hat) > 0 ) ) ) {
    print("All constraints were met!")
  }

  return(list(pred = f_X_hat, knots = knt_hat,
              coef = theta_hat, order = n) )
}

################################################################################
############################## CDF COEFFICIENTS ################################
################################################################################
theta_prime_func <- function(theta, knots, n) {
  p <- length(theta)

  if (length(knots) < p + n) {
    stop("`knots` must have length at least length(theta) + n.")
  }

  weights <- (knots[(n + 1):(p + n)] - knots[1:p]) / n

  c(0, cumsum(theta * weights))
}

validate_coef <- function(theta, B, y,
                          tol_end = 1e-6,
                          tol_fit = 1e-3) {

  p <- length(theta)

  # endpoints: CDF-style defaults theta[1]=0, theta[p]=1
  end_ok <- (abs(theta[1] - 0) <= tol_end) && (abs(theta[p] - 1) <= tol_end)

  # monotone
  d <- diff(theta)
  mono_ok <- all(d >= -1e-6)

  # fit error (RMS)
  pred <- as.numeric(B %*% theta)
  rms <- sqrt(mean((pred - y)^2))
  fit_ok <- is.finite(rms) && (rms <= tol_fit)

  list(ok = end_ok && mono_ok && fit_ok,
       end_ok = end_ok, mono_ok = mono_ok, fit_ok = fit_ok, rms = rms)
}

cdf_coef <- function(B, y,
                     tol_end = 1e-6,
                     tol_fit = 1e-3)
  {
  theta_fast <- cdf_coef_fast(B, y)

  # enforce endpoints for the fast attempt (cheap)
  theta_fast[1] <- 0
  theta_fast[length(theta_fast)] <- 1

  chk <- validate_coef(theta_fast, B, y, tol_end = tol_end, tol_fit = tol_fit)

  if (chk$ok) {
    return(theta_fast)
  }

  # fallback
  capture.output({
    theta <- cdf_coef_nloptr(B, y)$coef
  })
  return(theta)
}

# 1) Fast coefficient recovery
cdf_coef_fast <- function(B, y) {
  XtX <- crossprod(B)
  Xty <- crossprod(B, y)

  XtX_inv <- tryCatch(
    chol2inv(chol(XtX)),
    error = function(e1) tryCatch(
      solve(XtX),
      error = function(e2) MASS::ginv(XtX)
    )
  )
  as.numeric(XtX_inv %*% Xty)
}

# matcb <- crossprod(basisMatrix_F_X_hat)
# matcbinv <- tryCatch({
#   chol2inv(chol(matcb))  # Fastest if SPD
# }, error = function(e1) {
#   message("Matrix not SPD, using solve().")
#   tryCatch({
#     solve(matcb)
#   }, error = function(e2) {
#     message("Matrix singular, using ginv().")
#     MASS::ginv(matcb)
#   })
# })
# coef_F_X <- as.numeric(matcbinv %*% t(basisMatrix_F_X_hat) %*% F_X_hat)

# 2) Recover monotone cdf spline coefficients by solving a constrained nonlinear LS problem
# w.r.t to F_X_hat
cdf_coef_nloptr <- function(B, y) {

  p <- NCOL(B)

  # Define f_X, the density function
  F_X <- function(params) {
    theta <- params[1:p]
    return(B %*% theta)
  }

  # objFun <- function(params) {
  #   theta <- params[1:p]
  #
  #   # Compute F_values
  #   F_values <- pmax(B %*% theta, 1e-3)
  #
  #   # Log-likelihood computation
  #   rss <- sum((F_values - y)^2)
  #
  #   # Gradient w.r.t. theta
  #   grad_theta <- 2 * t(B) %*% (F_values - y)
  #
  #   return(list("objective" = rss, "gradient" = grad_theta))
  # }

  objFun <- function(params) {
    theta <- params[1:p]
    r <- as.numeric(B %*% theta - y)
    # Log-likelihood computation
    rss <- sum(r^2)
    # Gradient w.r.t. theta
    grad_theta <- 2 * as.numeric(crossprod(B, r))
    list(objective = rss, gradient = grad_theta)
  }

  # 2. Define constraints
  # 2.1. Inequality constraints
  ineq_constraints <- function(params) {
    # re-formulate ineq constraints to be of form g(x) <= 0
    theta <- params[1:p]

    # θ_1<...<θ_p
    epsilon <- 1e-6 # to ensure strict inequality
    constr2 <- c(epsilon - diff(theta))
    # Initialize the matrix with zeros
    grad2 <- matrix(0, nrow = p-1, ncol = p)
    for (i in 1:(p-1)) {
      grad2[i, i] <- 1
      grad2[i, i + 1] <- -1
    }


    return( list( "constraints" = constr2, "jacobian" = grad2 ) )
  }

  # 2.2. Equality constraints
  eq_constraints <- function(params) {
    theta <- params[1:p]
    # θ_1 = 0; θ_p = 0
    constr2 <- c(theta[1] - 0,
                 theta[p] - 1)
    # Initialize the matrix with zeros
    grad2 <- matrix(0, nrow = 2, ncol = length(params))
    # Set the first element of the first row to 1
    grad2[1, 1] <- 1
    # Set the last element of the last row to 1
    grad2[2, p] <- 1

    return( list( "constraints" = constr2, "jacobian" = grad2 ) )
  }

  # Initialize coefficients
  init_theta <- lm.fit(B, y)$coefficients
  init_theta[is.na(init_theta)] <- 1e-3
  init_theta <- pmin(pmax(init_theta, 1e-6), 1 - 1e-6)
  init_theta[1] <- 0
  init_theta[p] <- 1

  # Sort to enforce increasing trend before optimization starts
  init_theta <- sort(init_theta)

  # Perform Constrained Optimization using nloptr
  local_opts <- list(
    "algorithm" = "NLOPT_LD_MMA",
    "xtol_rel" = 1.0e-10
  )
  opts <- list(
    "algorithm" = "NLOPT_LD_AUGLAG",
    "xtol_rel" = 1.0e-10,
    "maxeval" = 10000,
    "local_opts" = local_opts
  )

  opt_result <- nloptr(
    x0 = init_theta,
    eval_f = objFun,
    eval_g_ineq = ineq_constraints,
    eval_g_eq = eq_constraints,
    opts = opts
  )

  # Extract results
  theta_hat <- opt_result$solution[1:p]

  y_prime <- B %*% theta_hat

  # Check
  if ( all(theta_hat >= 0) ) {
    print("All constraints were met!")
  }

  return(list(pred = y_prime, coef = theta_hat) )
}

theta_prime_bivariate_func <- function(theta, knots1, knots2, n1, n2) {
  p1 <- length(knots1) - n1
  p2 <- length(knots2) - n2

  if (length(theta) != p1 * p2) {
    stop("length(theta) must be equal to (length(knots1)-n1) * (length(knots2)-n2).")
  }

  # reshape theta vector into matrix
  theta_mat <- matrix(theta, nrow = p1, ncol = p2, byrow = TRUE)

  # weights
  w1 <- (knots1[(1 + n1):(p1 + n1)] - knots1[1:p1]) / n1
  w2 <- (knots2[(1 + n2):(p2 + n2)] - knots2[1:p2]) / n2

  # elementwise weighting
  weighted_theta <- theta_mat * outer(w1, w2)

  # cumulative sums over both dimensions
  theta_prime <- apply(weighted_theta, 2, cumsum)
  theta_prime <- t(apply(theta_prime, 1, cumsum))

  # add zero row and zero column
  theta_prime <- rbind(0, theta_prime)
  theta_prime <- cbind(0, theta_prime)

  as.vector(t(theta_prime))
}


# Convert tensor-product density coefficients to joint-CDF coefficients.
# Coefficients and knots follow GeDS tensor ordering (last dimension fastest).
theta_prime_multivariate <- function(theta, knots, orders) {
  if (!is.list(knots) || !length(knots)) {
    stop("'knots' must be a non-empty list.", call. = FALSE)
  }
  if (length(orders) == 1L) orders <- rep(orders, length(knots))
  if (length(orders) != length(knots)) {
    stop("'orders' must have one value per dimension.", call. = FALSE)
  }

  basis_counts <- lengths(knots) - orders
  if (any(basis_counts < 1L) || length(theta) != prod(basis_counts)) {
    stop("The coefficient count does not match the supplied knots and orders.",
         call. = FALSE)
  }
  weights <- Map(function(k, n) diff(k, lag = n) / n, knots, orders)
  weighted <- as.numeric(theta) * ddfs_tensor_weights(weights)

  # Reversed dimensions map R column-major arrays to GeDS tensor ordering.
  cumulative <- array(weighted, dim = rev(basis_counts))
  for (margin in seq_along(dim(cumulative))) {
    permutation <- c(margin, setdiff(seq_along(dim(cumulative)), margin))
    inverse <- order(permutation)
    permuted <- aperm(cumulative, permutation)
    dims <- dim(permuted)
    permuted <- apply(matrix(permuted, nrow = dims[1L]), 2L, cumsum)
    dim(permuted) <- dims
    cumulative <- aperm(permuted, inverse)
  }

  padded <- array(0, dim = dim(cumulative) + 1L)
  target <- lapply(dim(padded), function(size) 2L:size)
  padded <- do.call(`[<-`, c(list(padded), target, list(value = cumulative)))
  result <- as.vector(padded)
  attr(result, "basis.index") <- ddfs_tensor_index(
    basis_counts + 1L, names(knots)
  )
  result
}
