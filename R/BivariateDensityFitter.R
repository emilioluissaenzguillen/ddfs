################################################################################
################################################################################
########################### Bivariate Density Fitter ###########################
################################################################################
################################################################################
#' @importFrom MASS ginv
#' @importFrom graphics par legend
#' @importFrom plot3D scatter3D
#' @importFrom splines splineDesign
#' @importFrom utils getFromNamespace
BivariateDensityFitter <- function(XY, n = 3L, min_iterations = 1,
                                   max_iterations = 50L, max.intknots = 1,
                                   beta = 0.1, phi_F_XY = 0.95, q_F_XY = 2,
                                   tails_count_threshold = 0.035, plot = FALSE, resids_plot = FALSE,
                                   pdf = NULL, cdf = NULL, stop_plotting = 0)
{
  n <- as.integer(n)
  N <- NROW(XY)

  # Initialize
  X <- XY[,1]; Y <- XY[,2]
  F_XY <- apply(XY, 1, function(v) bivariate_ecdf(XY, v[1], v[2]))
  resid_XY <- F_XY
  Xextr <- range(X); Yextr <- range(Y)
  X_grid <- seq(from = min(X), to = max(X), length.out = round(sqrt(N)))
  Y_grid <- seq(from = min(Y), to = max(Y), length.out = round(sqrt(N)))
  grid_data <- expand.grid(X = X_grid, Y = Y_grid)

  f_XY_hat <- F_XY_hat <- InterKnotsX <- InterKnotsY <-  NULL
  f_XY_hat_list <- F_XY_hat_list <- list()
  RSS <- numeric(0)

  out <- list(
    f_XY_hat = f_XY_hat,
    F_XY_hat = F_XY_hat, Type = "Biv - DDFS",
    Args = list(XY = XY, ecdf = F_XY, phi = phi_F_XY, q = q_F_XY, beta = beta),
    model = NULL
    )

  # Calculate proportion of values in the lower and upper 5% of the range
  count_lowerX <- sum(X <= range(X)[1] + diff(range(X)) * 0.05) / length(X) # lower 5%
  count_upperX <- sum(X > range(X)[2] - diff(range(X)) * 0.05) / length(X)  # upper 5%
  count_lowerY <- sum(Y <= range(Y)[1] + diff(range(Y)) * 0.05) / length(Y) # lower 5%
  count_upperY <- sum(Y > range(Y)[2] - diff(range(Y)) * 0.05) / length(Y)  # upper 5%
  threshold <- tails_count_threshold

  leftX_decreasingtail <- count_lowerX < threshold
  rightX_decreasingtail <- count_upperX < threshold
  leftY_decreasingtail <- count_lowerY < threshold
  rightY_decreasingtail <- count_upperY < threshold

  # print(paste0("leftX_decreasingtail=", leftX_decreasingtail, " ", count_lowerX,
  #              " rightX_decreasingtail=", rightX_decreasingtail, " ", count_upperX))
  # print(paste0("leftY_decreasingtail=", leftY_decreasingtail, " ", count_lowerY,
  #              " rightY_decreasingtail=", rightY_decreasingtail, " ", count_upperY))

  # Iterate
  for (iter in 1:max_iterations) {

    # 2) Fit a linear GeDS regression model with one internal knot to the ECDF of the data.
    if (iter > 1) {
      suppressWarnings({
        Gmod <- NGeDS(resid_XY ~ f(X,Y), beta = beta,
                      max.intknots = max.intknots + iter - 2,
                      intknots_init = list(ikX = InterKnotsX, ikY = InterKnotsY),
                      higher_order = FALSE)
      })

      if ( length(Gmod$linear.intknots$Xk) == length(InterKnotsX) &&  length(Gmod$linear.intknots$Yk) == length(InterKnotsY) ) {
        break("No more knots being added.")
      }
      InterKnotsX <- Gmod$linear.intknots$Xk
      InterKnotsY <- Gmod$linear.intknots$Yk
    }

    kntX <- sort(c(InterKnotsX,rep(Xextr,n)))
    kntY <- sort(c(InterKnotsY,rep(Yextr,n)))
    p1 <- length(InterKnotsX) + n; p2 <- length(InterKnotsY) + n
    p <- p1 * p2

    ##################
    # Residuals plot #
    ##################
    if(plot && iter > 1 && resids_plot) {

      par(mai = c(0.42, 0, 0.52, 0)) # bottom, left, top, right par(mai = c(0.52, 0, 0.52, 0)) # par(mai = c(1.02, 0.82, 0.82, 0.42))
      plot(Gmod, n = 2, xlab = "", ylab = "", zlab = "", main = "detail",
           legend.text = c(expression(rho >= hat(rho)), expression(rho < hat(rho))) )

      text3D(
        range(X)[1] - 0.2 * diff(range(X)),
        range(Y)[1] - 0.2 * diff(range(Y)),
        sum(range(predict(Gmod, n = 2, newdata = grid_data)))/2,
        labels = expression(rho),
        add = TRUE,
        cex = 1
      )
      text3D(round(sum(range(X))/2), range(Y)[1] - 0.15*diff(range(Y)), min(predict(Gmod, n = 2, newdata = grid_data)),
             expression(X[1]), add = TRUE, cex = 1, srt = -45)
      text3D(range(X)[2] + 0.15*diff(range(X)), round(sum(range(Y))/2), min(predict(Gmod, n = 2, newdata = grid_data)),
             expression(X[2]), add = TRUE, cex = 1, srt = 45)
      par(mai = c(1.02, 0.82, 0.82, 0.42))
    }

    # Helper to produce plots
    if (stop_plotting !=0 && stop_plotting + 1 == iter) break

    # 3) Constrained MLE estimation of coefficients
    # 3.1) Initialize coefficients
    # \theta_j^0=\frac{n}{\sum_{j=1}^p(t_j-t_{j-n})}

    # Conditionally initialize theta based on the tail flags
    if (leftX_decreasingtail) leftX = 1 else leftX = NULL
    if (rightX_decreasingtail) rightX = p1 else rightX = NULL
    if (leftY_decreasingtail) leftY = 1 else leftY = NULL
    if (rightY_decreasingtail) rightY = p2 else rightY = NULL

    if (is.null(leftX) && is.null(rightX)) diffX <- diff(kntX, n) else diffX <- diff(kntX, n)[-c(leftX, rightX)]
    if (is.null(leftY) && is.null(rightY)) diffY <- diff(kntY, n) else diffY <- diff(kntY, n)[-c(leftY, rightY)]

    theta <- rep(n^2/sum(outer(diffX, diffY, `*`)), p)
    theta_matrix <- matrix(theta, nrow = p1, ncol = p2)

    # Conditionally set boundaries to zero based on the tail flags
    if (leftX_decreasingtail)   theta_matrix[1, ]             <- 0  # Lower X
    if (rightX_decreasingtail)  theta_matrix[nrow(theta_matrix), ] <- 0  # Upper X
    if (leftY_decreasingtail)   theta_matrix[, 1]             <- 0  # Lower Y
    if (rightY_decreasingtail)  theta_matrix[, ncol(theta_matrix)] <- 0  # Upper Y

    # Flatten back to vector form
    theta <- as.vector(t(theta_matrix))
    theta_prev <- theta
    indexes <- which(theta != 0)

    N_init <- N - sum(c(leftX_decreasingtail, rightX_decreasingtail, leftY_decreasingtail, rightY_decreasingtail))
    first_term <- as.vector(t( (1/N_init) * (n^2 / outer(diff(kntX, n), diff(kntY, n), `*`)) ))


    # Create spline basis matrix using specified knots, evaluation points and order
    basisMatrixX <- splineDesign(knots = kntX, x = X, ord = n,
                                 derivs = rep(0,length(X)), outer.ok = T)
    basisMatrixY <- splineDesign(knots = kntY, x = Y, ord = n,
                                 derivs = rep(0,length(Y)), outer.ok = T)
    basisMatrixbiv <- getFromNamespace("tensorProd", "GeDS")(basisMatrixX, basisMatrixY)


    theta <- constrMLE_R(theta, indexes, basisMatrixbiv, first_term)

    # thetaarma <- as.numeric(constrMLE_arma(theta, indexes, basisMatrixbiv, first_term))
    #
    # if ( all(abs(theta - thetaarma) < 1e-5) ) print("Both equal!")
    #
    # library(microbenchmark)
    # # Run the benchmark
    # mb_results <- microbenchmark(
    #   R_version = constrMLE_R(theta, indexes, basisMatrixbiv, first_term),
    #   arma_version = constrMLE_arma(theta, indexes, basisMatrixbiv, first_term),
    #   times = 1000
    # )
    # # Extract the median execution times
    # medians <- summary(mb_results)$median
    # names(medians) <- summary(mb_results)$expr
    # # Find the faster version
    # faster_version <- names(medians)[which.min(medians)]
    # # Print results
    # cat("Median execution times (in nanoseconds):\n")
    # print(medians)
    # cat("\nThe faster implementation is:", faster_version, "\n")


    # Calculate the difference
    diff <- abs( sum(theta / (N_init * first_term) ) - 1)
    if(diff !=0) {
      # Calculate the adjustment factor
      adjustment_factor <- 1 - diff
      # Adjust theta uniformly
      theta <- theta / adjustment_factor
    }

    # Check
    diff <- abs(sum(theta * as.vector(t(outer(diff(kntX, n), diff(kntY, n), `*`))) / n^2) - 1)
    if( !all( all(theta >= 0) && diff < 1e-6 ) ) {
      warning(paste0("Constraints were not met by ", diff))
    }

    f_XY_hat <- basisMatrixbiv %*% theta

    # Save the current f_XY_hat_list
    f_XY_hat_list[[iter]] <- list(pred = f_XY_hat, knots = list(Xk = kntX, Yk = kntY),
                                  coef = theta, order = n)

    ################
    # Density plot #
    ################
    if(plot) {
      out$f_XY_hat <- f_XY_hat_list[[iter]]
      par(mai = c(0.42, 0, 0.52, 0)) # bottom, left, top, right
      plot.ddfs(out, type = "density", f = pdf)
      par(mai = c(1.02, 0.82, 0.82, 0.42))
    }

    # 4) Based on the latter estimated density, we obtain a corresponding estimate of the CDF
    # c.f. Dierckx (1993), Chapter 2, formula (20)
    # Define BivariateSpline object
    F_XY_hat <- compute_bivariate_integral(X = X, Y = Y, kntX = kntX, kntY = kntY,
                                               theta = theta, n = n)

    resid_XY <- F_XY - F_XY_hat
    RSS[iter] <- sum(resid_XY^2)

    # To calculate coefficients much better to use grid_data
    F_XY_hat_grid <- compute_bivariate_integral(X = grid_data[,1], Y = grid_data[,2],
                                                kntX = kntX, kntY = kntY, theta = theta, n = n)

    # The integral of a spline of degree n is a spline of degree n + 1
    basisMatrix_F_XY_hatX <- splineDesign(knots = sort(c(InterKnotsX,rep(Xextr,n+1))),
                                          x = grid_data[,1], ord = n+1, derivs = rep(0,length(grid_data[,1])),
                                          outer.ok = T)
    basisMatrix_F_XY_hatY <- splineDesign(knots = sort(c(InterKnotsY,rep(Yextr,n+1))),
                                          x = grid_data[,2], ord = n+1, derivs = rep(0,length(grid_data[,2])),
                                          outer.ok = T)
    basisMatrixbiv_F_XY_hat <- getFromNamespace("tensorProd", "GeDS")(basisMatrix_F_XY_hatX, basisMatrix_F_XY_hatY)

    # Coefficients of F_XY_hat
    matcb <- crossprod(basisMatrixbiv_F_XY_hat)
    matcbinv <- tryCatch({
      chol2inv(chol(matcb))  # Fastest if SPD
    }, error = function(e1) {
      message("Matrix not SPD, using solve().")
      tryCatch({
        solve(matcb)
      }, error = function(e2) {
        message("Matrix singular, using ginv().")
        MASS::ginv(matcb)
      })
    })

    # Save the current F_XY_hat_list
    F_XY_hat_list[[iter]] <- list(pred = F_XY_hat,
                                  knots = list(Xk = sort(c(InterKnotsX,rep(Xextr,n+1))),
                                               Yk = sort(c(InterKnotsY,rep(Yextr,n+1)))),
                                  coef = matcbinv %*% t(basisMatrixbiv_F_XY_hat) %*% F_XY_hat_grid,
                                  order = n+1)

    # print(matrix(F_XY_hat_list[[iter]]$coef, nrow = p1+1, ncol = p2+1, byrow=FALSE))
    # print(as.numeric(round(basisMatrixbiv_F_XY_hat %*% F_XY_hat_list[[iter]]$coef - F_XY_hat_grid, 4)))

    #################
    # Integral plot #
    #################
    if(plot) {
      out$F_XY_hat <- F_XY_hat_list[[iter]]
      par(mai = c(0.42, 0, 0.52, 0)) # bottom, left, top, right
      plot.ddfs(out, type = "distribution", f = cdf)
      par(mai = c(1.02, 0.82, 0.82, 0.42))
    }

    #################
    # Stopping rule #
    #################
    if (iter > q_F_XY && iter > min_iterations) {
      if (RSS[iter]/RSS[iter - q_F_XY] >= phi_F_XY) {
        #cat("Stopping iterations due to small improvement in deviance\n\n")
        break
      }
    }

  }

  ######################
  ### Residuals Plot ###
  ######################
  if (plot && resids_plot && stop_plotting == max_iterations && stop_plotting == iter) {
    # Plot the perspective 3D surface defined by X_grid, Y_grid, and the fitted values
    par(mai = c(0.42, 0, 0.52, 0))
    scatter3D(x = X, y = Y, z = resid_XY, phi = 25, theta = 50,
              xlab = '', ylab = '', zlab = "residuals", main = "", zlim = range(resid_XY),
              ticktype = "detailed", expand = 0.5, colkey = FALSE, border = "black")
    # Add title & labels
    text3D(0.25, -4, range(resid_XY)[1], expression(X[1]), add = TRUE, cex = 1, srt = -45)
    text3D(4, -0.4, range(resid_XY)[1], expression(X[2]), add = TRUE, cex = 1, srt = 45)

    legend("topright",
           legend = c(expression(rho == F[N](x) - hat(F)(x))),
           col = "black",
           pch = 19,
           bty = "n")
    # Reset to default margins and mgp
    par(mai = c(0.42, 0, 0.52, 0))
  }


  if (iter != max_iterations) {
    index <- iter-q_F_XY
  } else {
    index <- iter
  }

  out$f_XY_hat <- f_XY_hat_list[[index]]
  out$F_XY_hat <- F_XY_hat_list[[index]]
  out$RSS$mindist <- RSS
  out$model <- index

  return(out)
}



