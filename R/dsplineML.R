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

cdf_coef <- function(basisMatrix_F_X_hat, F_X_hat) {


  p <- NCOL(basisMatrix_F_X_hat)

  # Define f_X, the density function
  F_X <- function(params) {
    theta <- params[1:p]
    return(basisMatrix_F_X_hat %*% theta)
  }

  objFun <- function(params) {
    theta <- params[1:p]

    # Compute F_values
    F_values <- pmax(basisMatrix_F_X_hat %*% theta, 1e-3)

    # Log-likelihood computation
    rss <- sum((F_values - F_X_hat)^2)

    # Gradient w.r.t. theta
    grad_theta <- 2 * t(basisMatrix_F_X_hat) %*% (F_values - F_X_hat)

    return(list("objective" = rss, "gradient" = grad_theta))
  }

  # 2. Define constraints
  # 2.1. Inequality constraints
  ineq_constraints <- function(params) {
    # re-formulate ineq constraints to be of form g(x) <= 0
    theta <- params[1:p]

    # (2) θ_1<...<θ_p
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
    # (2) θ_1 = 0; θ_p = 0
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
  init_theta <- lm.fit(basisMatrix_F_X_hat, F_X_hat)$coefficients
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

  F_X_hat_prime <- basisMatrix_F_X_hat %*% theta_hat

  # Check
  if ( all(theta_hat >= 0) ) {
    print("All constraints were met!")
  }

  return(list(pred = F_X_hat_prime, coef = theta_hat) )
}

