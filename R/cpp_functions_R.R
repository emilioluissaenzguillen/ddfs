constrMLE_R <- function(theta, indexes, basisMatrix, first_term, tol = 1e-6) {
  # Initialize theta_prev as a copy of theta
  theta_prev <- theta

  ok <- TRUE
  while (ok) {
    theta_prev[indexes] <- theta[indexes]
    f_hat_prev <- basisMatrix %*% theta_prev
    f_hat_prev[f_hat_prev == 0] <- 1e-6
    theta <- first_term * theta_prev * colSums(basisMatrix / as.vector(f_hat_prev))

    # Convergence check: if the maximum change for the indexes is less than 1e-6, exit the loop
    if (max(abs(theta[indexes] - theta_prev[indexes])) < tol) ok <- FALSE
  }

  return(theta)
}


