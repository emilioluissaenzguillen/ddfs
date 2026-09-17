################################################################################
################################################################################
########################## Univariate Density Fitter ###########################
################################################################################
################################################################################
#' @importFrom diptest dip.test
#' @importFrom GeDS NGeDS PPolyRep
#' @importFrom graphics hist
#' @importFrom splines splineDesign
#' @importFrom splines2 ibs
#' @importFrom utils capture.output

UnivariateDensityFitter <- function(X, n = 4L, min_iterations = 2,
                                    max_iterations = 50L, max.intknots = 1,
                                    beta = 0, phi_F_X = 0.3, q_F_X = 1, stoptype = "RDMD",
                                    tail_decay = c("auto", "none", "left", "right", "both"),
                                    boundary_extension = 0,
                                    plot = FALSE, resids_plot = FALSE, pdf = NULL, cdf = NULL,
                                    stop_plotting = 0, schoenberg = FALSE)
{

  n <- as.integer(n)
  if (anyNA(X)) {
    warning("Missing values (NA) in 'X' have been removed.")
    X <- X[!is.na(X)]
  }
  X <- sort(as.numeric(X))
  N <- NROW(X)

  # Initialize
  F_X <- cumsum(rep(1/N, N))
  if (length(boundary_extension) == 1L) {
    boundary_extension <- rep(boundary_extension, 2L)
  }
  extension <- boundary_extension * stats::bw.nrd0(X)
  extr <- range(X) + c(-extension[1L], extension[2L])
  f_X_hat <- F_X_hat <- InterKnots <- NULL
  oldintc <- oldslp <- phis <- phis_star <- NULL
  f_X_hat_list <- F_X_hat_list <- list()
  RSS <- numeric(0); RSS_GeDS <- numeric(0)

  out <- list(
    f_X_hat = f_X_hat,
    F_X_hat = F_X_hat, type = "Univ - DDFS",
    args = list(X = X, ecdf = F_X, phi = phi_F_X, q = q_F_X, beta = beta,
                boundary_extension = boundary_extension),
    RSS = list(mindist = NULL, GeDS = NULL),
    model = NULL
  )

  tail_decay <- match.arg(tail_decay)

  if (tail_decay == "auto"){
    # # Option 1
    # t_function <- function(N) log(N)/N
    #
    # count_lower <- sum(X <= range(X)[1] + diff(range(X)) * t_function(length(X)) ) / length(X) # lower 5%; length(X)/10000
    # count_upper <- sum(X > range(X)[2] - diff(range(X)) * t_function(length(X)) ) / length(X)  # upper 5%
    # threshold <- tails_count_threshold
    #
    # left_decreasingtail <- count_lower < threshold
    # right_decreasingtail <- count_upper < threshold


    # Option 2
    counts <- hist(X, plot = FALSE)$counts # breaks = sqrt(length(X))*(1+length(X)^(-0.1))

    if (suppressMessages(suppressWarnings(dip.test(X)$p.value)) >= 0.3) {
      Gmod <- suppressMessages(suppressWarnings(NGeDS(counts ~ f(seq_along(counts)),
                                                      min.intknots = 1, max.intknots = 1)))
    } else {
      Gmod <- suppressMessages(suppressWarnings(NGeDS(counts ~ f(seq_along(counts)),
                                                      phi = 0.99)))
    }

    # plot(Gmod, n = 2)
    Gmod <- PPolyRep(Gmod, n = 2)

    left_decreasingtail <- right_decreasingtail <- TRUE

    if (Gmod$coefficients[1,2] < 5/2) {
      left_decreasingtail <- FALSE
    }

    if (Gmod$coefficients[nrow(Gmod$coefficients)-1, 2] > -5/2) {
      right_decreasingtail <- FALSE
    }


    # # Option 3
    # # Check whether f(X) → 0 as X → ±∞
    # # Calculate proportion of values in the lower and upper 5% of the range
    # t_function <- function(N) {
    #   exp(1/3)*log(N) / (N^(1 + 1/log(N)))
    # }
    # count_lower <- sum(X <= range(X)[1] + diff(range(X)) * t_function(length(X)) ) / length(X) # lower 5%; length(X)/10000
    # count_upper <- sum(X > range(X)[2] - diff(range(X)) * t_function(length(X)) ) / length(X)  # upper 5%
    # threshold <- tails_count_threshold * (1+exp(1)/length(X)^(1/3))
    #
    # left_decreasingtail <- count_lower < threshold
    # right_decreasingtail <- count_upper < threshold
    #
    # # Extra tail check
    # counts <- hist(X, plot = FALSE)$counts # breaks = sqrt(length(X))*(1+length(X)^(-0.1))
    #
    # if (suppressMessages(suppressWarnings(dip.test(X)$p.value)) >= 0.3) {
    #   Gmod <- suppressMessages(suppressWarnings(NGeDS(counts ~ f(seq_along(counts)),
    #                                                   min.intknots = 1, max.intknots = 1)))
    # } else {
    #   Gmod <- suppressMessages(suppressWarnings(NGeDS(counts ~ f(seq_along(counts)),
    #                                                   phi = 0.99)))
    # }
    #
    # # plot(Gmod, n = 2)
    # Gmod <- PPolyRep(Gmod, n = 2)
    #
    # if (Gmod$coefficients[1,2] < 1/3) {
    #   left_decreasingtail <- FALSE
    # } else if (Gmod$coefficients[nrow(Gmod$coefficients)-1, 2] > -1/3) {
    #   right_decreasingtail <- FALSE
    # }

    # # Check whether f(X) → 0 as X → ±∞
    # # Calculate proportion of values in the lower and upper 5% of the range
    # count_lower <- sum(X <= range(X)[1] + diff(range(X)) * 0.05) / length(X) # lower 5%
    # count_upper <- sum(X > range(X)[2] - diff(range(X)) * 0.05) / length(X)  # upper 5%
    # threshold <- tails_count_threshold
    #
    # left_decreasingtail <- count_lower < threshold
    # right_decreasingtail <- count_upper < threshold

    # left_decreasingtail <- TRUE; right_decreasingtail <- TRUE

    # print(paste0("left_decreasingtail=", left_decreasingtail, " right_decreasingtail=", right_decreasingtail))

  } else if (tail_decay == "both") {
    left_decreasingtail <- right_decreasingtail <- TRUE
  } else if (tail_decay == "none") {
    left_decreasingtail <- right_decreasingtail <- FALSE
  } else if (tail_decay == "left") {
    left_decreasingtail <- TRUE; right_decreasingtail <- FALSE
  } else if (tail_decay == "right") {
    left_decreasingtail <- FALSE; right_decreasingtail <- TRUE
  }




  # Iterate
  for (iter in 1:max_iterations) {

    # 2) Fit a linear GeDS regression model with one internal knot to the ECDF of the data.
    if (iter > 1) {

      if (!schoenberg) {
      Gmod <- suppressMessages(suppressWarnings(
        NGeDS(resid_X ~ f(X), beta = beta, max.intknots = max.intknots + iter - 2,
                      intknots_init = InterKnots, higher_order = FALSE)
        ))
      X <- Gmod$args$X

      if ((is.null(Gmod$linear.intknots)  || length(Gmod$linear.intknots) == length(InterKnots)) && iter > 1) {
        iter <- iter - 1; max_iterations <- iter; break("No more knots being added.")
        }
      InterKnots <- Gmod$linear.intknots

      } else {
        Gmod1 <- suppressMessages(suppressWarnings(
          NGeDS(resid_X ~ f(X), beta = beta, max.intknots = max.intknots + iter - 2,
                intknots_init = InterKnots, higher_order = FALSE)
        ))
        Gmod2 <- suppressMessages(suppressWarnings(
          NGeDS(resid_X ~ f(X), beta = beta, max.intknots = max.intknots + iter - 1,
                intknots_init = InterKnots, higher_order = TRUE)
        ))
        Gmod3 <- suppressMessages(suppressWarnings(
          NGeDS(resid_X ~ f(X), beta = beta, max.intknots = max.intknots + iter,
                intknots_init = InterKnots, higher_order = TRUE)
        ))

        # Create a named vector of deviances
        devs <- c(linear = Gmod1$dev.linear,
                  quadratic = Gmod2$dev.quadratic,
                  cubic = Gmod3$dev.cubic)
        # Find model with the lowest deviance
        best_model <- names(which.min(devs))
        # cat("The model with the lowest deviance is:", best_model, "with deviance =", devs[best_model], "\n")

        if (best_model == "linear") {
          X <- Gmod1$args$X
          if (is.null(Gmod1$linear.intknots) && iter > 1) {
            iter <- iter - 1; max_iterations <- iter; break("No more knots being added.")
            }
          InterKnots <- Gmod1$linear.intknots
          Gmod <- Gmod1

          } else if (best_model == "quadratic") {
            X <- Gmod2$args$X
            if (is.null(Gmod2$quadratic.intknots) && iter > 1) {
              iter <- iter - 1; max_iterations <- iter; break("No more knots being added.")
              }
            InterKnots <- Gmod2$quadratic.intknots
            Gmod <- Gmod2

          } else if (best_model == "cubic") {
            X <- Gmod3$args$X
            if (is.null(Gmod3$cubic.intknots) && iter > 1) {
              iter <- iter - 1; max_iterations <- iter; break("No more knots being added.")
            }
            InterKnots <- Gmod3$cubic.intknots
            Gmod <- Gmod3

          }

      }

    }

    ##################
    # Residuals plot #
    ##################
    if(plot && iter > 1 && resids_plot) {
      # To keep the ylim constant
      if (iter == 2) resid_X_1 <- resid_X
      if (!schoenberg) {
        ord = 2
      } else {
        if (best_model == "linear") {
          ord = 2
          } else if (best_model == "quadratic") {
            ord = 3
            } else if (best_model == "cubic") {
              ord = 4
            }

      }
      suppressWarnings(
        plot(Gmod, n = ord, ylim = range(resid_X_1, Gmod$linear.fit$predicted),
             xlab = "", ylab = "residuals",
             main = "detail", legend.pos = "topleft",
             legend.text = c(expression(rho == F[N](x) - hat(F)[DDFS](x)), "GeDS"))
      )

    }

    # Helper to produce plots
    if (stop_plotting !=0 && stop_plotting + 1 == iter) break

    # 3) Constrained MLE estimation of coefficients
    knt <- sort(c(InterKnots,rep(extr,n)))
    p <- length(InterKnots) + n


    # Conditional logic based on the density of points in the lower and upper 5%
    if ( !left_decreasingtail && !right_decreasingtail ) {
      # Case: High density in both tails
      # 3.1) Initialize coefficients
      # \theta_j^0=\frac{n}{\sum_{j=1}^p(t_j-t_{j-n})}
      theta <- rep(n / sum(diff(knt, n)), p)
      theta_prev <- theta
      # 3.2)
      # sum_j <- 0
      # for (j in (n+1):(p+n)) {
      #   sum_j <- sum_j + (knt[j] - knt[j - n])
      # }
      first_term <- (1 / N) * (n / diff(knt, n))
      indexes <- 1:p

    } else if ( left_decreasingtail && right_decreasingtail ) {
      # Case: Low density in both tails
      theta <- c(0, rep(n / sum(diff(knt, n)[-c(1, p)]), p - 2), 0)
      theta_prev <- theta
      boundary_count <- sum(X == extr[1L]) + sum(X == extr[2L])
      first_term <- (1 / (N - boundary_count)) * (n / diff(knt, n))
      indexes <- 2:(p - 1)

    } else if ( !left_decreasingtail && right_decreasingtail ) {
      # Case: Low density in the upper tail only
      theta <- c(rep(n / sum(diff(knt, n)[-p]), p - 1), 0)
      theta_prev <- theta
      first_term <- (1 / (N - sum(X == extr[2L]))) * (n / diff(knt, n))
      indexes <- 1:(p - 1)

    } else if ( left_decreasingtail && !right_decreasingtail ) {
      # Case: Low density in the lower tail only
      theta <- c(0, rep(n / sum(diff(knt, n)[-1]), p - 1))
      theta_prev <- theta
      first_term <- (1 / (N - sum(X == extr[1L]))) * (n / diff(knt, n))
      indexes <- 2:p
    }

    # Create spline basis matrix using specified knots, evaluation points and order
    basisMatrix <- splineDesign(knots = knt, x = X, ord = n, derivs = rep(0,length(X)),
                                outer.ok = T)

    # thetacpp <- constrMLE(theta, indexes, basisMatrix, first_term)
    # R function is the fastest one

    # thetaarma <- constrMLE_arma(theta, indexes, basisMatrix, first_term)
    if (N > 10^4) tol = 1e-3 else tol = 1e-6
    theta <- constrMLE_R(theta, indexes, basisMatrix, first_term, tol)

    # if ( all(abs(theta - thetaarma) < 1e-5) ) print("Both equal!")

    # library(microbenchmark)
    # # Run the benchmark
    # mb_results <- microbenchmark(
    #   R_version = constrMLE_R(theta, indexes, basisMatrix, first_term),
    #   arma_version = constrMLE_arma(theta, indexes, basisMatrix, first_term),
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


    # Check
    if( !all( all(theta >= 0) && abs(sum(theta * diff(knt, n)/n) - 1) < 1e-6 ) ) warning("Constraints were not met.")

    f_X_hat <- basisMatrix %*% theta

    # Save the current f_X_hat_list
    f_X_hat_list[[iter]] <- list(pred = f_X_hat, knots = knt, coef = theta, order = n)

    ################
    # Density plot #
    ################
    if (plot) {
      out$f_X_hat <- f_X_hat_list[[iter]]
      plot.ddfs(out, fit = "pdf", f = pdf)
    }

    # 4) Based on the latter estimated density, we obtain a corresponding estimate of the CDF
    # c.f. De Boor, 2001, Chapter X, formula (33); the integral of a spline of degree n is a spline of degree n + 1

    # F_X_hat <- Integrate(knots = knt, coef = theta, from = knt[1], to = X, n = n)
    # F_X_hat_ibs <- ibs::ibs(X, knt, ord = n, coef = theta)
    # F_X_hat <- ibs(x = X, degree = n - 1,
    #                knots = InterKnots, Boundary.knots = extr,
    #                coef = theta, intercept = TRUE) %*% theta

    basisMatrix_F_X_hat <- splineDesign(knots = sort(c(InterKnots,rep(extr,n+1))),
                                        x = X, ord = n+1, derivs = rep(0,length(X)),
                                        outer.ok = T)
    theta_prime <- theta_prime_func(theta, sort(c(InterKnots,rep(extr,n))), n)
    # coef_F_X <- cdf_coef(basisMatrix_F_X_hat, F_X_hat)
    # print(round(as.numeric(theta_prime - coef_F_X),5))
    F_X_hat <- basisMatrix_F_X_hat%*% theta_prime

    resid_X <- F_X - F_X_hat
    RSS[iter] <- sum(resid_X^2)
    if(iter > 1) {
      if (iter == 2) RSS_GeDS[1] <- Gmod$rss[1]
      RSS_GeDS[iter] <- Gmod$linear.fit$rss
    }

    # print(as.numeric(round(F_X_hat - basisMatrix_F_X_hat%*% theta_prime,4)))

    # Save the current F_X_hat_list
    F_X_hat_list[[iter]] <- list(pred = F_X_hat, knots = sort(c(InterKnots,rep(extr,n+1))),
                                 coef = theta_prime, order = n+1)

    #################
    # Integral plot #
    #################
    if(plot) {
      out$F_X_hat <- F_X_hat_list[[iter]]
      plot.ddfs(out, fit = "cdf", f = cdf)
    }

    ###################
    ## Stopping rule ##
    ###################
    ## 1) Minimum distance F_N - F_hat
    if ((stoptype == "RDMD"|| stoptype == "SRMD") && iter > q_F_X) {

      phi_temp <- RSS[iter]/RSS[iter - q_F_X]
      if (phi_temp >= 1) phi_temp <- 1-1e-6
      phis <- c(phis, phi_temp)

      if(iter > min_iterations) {
        # RD-MinDist stopping rule
        if (stoptype == "RDMD" && phi_temp >= phi_F_X) {
          #cat("Stopping iterations due to small improvement in RDMD\n\n")
          break
          # SR-MinDist stopping rule
        } else if (stoptype == "SRMD") {
          # \hat{φ}_κ = 1 − exp{\hat{γ}_0 + \hat{γ}_1*κ}
          # 1-\hat{φ}_κ = exp{\hat{γ}_0 + \hat{γ}_1*κ}
          # ln(1-\hat{φ}_κ) = \hat{γ}_0 + \hat{γ}_1*κ
          # Fit a linear model ln(1-φ) ~ \hat{γ}_0 + \hat{γ}_1*κ to the sample {φ_h, h}^κ_{h=q}
          phismod <- log(1-phis); kappa <- length(InterKnots)
          gamma <- .lm.fit(cbind(1, q_F_X:kappa), phismod)$coef
          # Calculate \hat{φ}_κ based on the estimated coefficients
          phi_kappa <- 1 - exp(gamma[1])*exp(gamma[2]*kappa)
          # Store \hat{φ}_κ and the estimated coefficients \hat{γ}_0 and \hat{γ}_1
          phis_star <- c(phis_star, phi_kappa)
          oldintc <- c(oldintc, gamma[1]); oldslp <- c(oldslp, gamma[2])
          # Check if \hat{φ}_κ ≥ φ_{exit}
          if (phi_kappa >= phi_F_X) {
            #cat("Stopping iterations due to small improvement in SRMD\n\n")
            break
          }
        }
      }

      ## 2) GeDS stopping rule (w.r.t. r_i)
    } else if ( (stoptype == "RD"|| stoptype == "SR") && iter > q_F_X) {

      phi_temp <- RSS_GeDS[iter]/RSS_GeDS[iter - q_F_X]
      if (phi_temp >= 1) phi_temp <- 1-1e-6
      phis <- c(phis, phi_temp)

      if(iter > min_iterations) {
        # GeDS-RD stopping rule
        if (stoptype == "RD" && phi_temp >= phi_F_X) {
          #cat("Stopping iterations due to small improvement in RD\n\n")
          break

          # GeDS-SR stopping rule
        } else if (stoptype == "SR") {
          # \hat{φ}_κ = 1 − exp{\hat{γ}_0 + \hat{γ}_1*κ}
          # 1-\hat{φ}_κ = exp{\hat{γ}_0 + \hat{γ}_1*κ}
          # ln(1-\hat{φ}_κ) = \hat{γ}_0 + \hat{γ}_1*κ
          # Fit a linear model ln(1-φ) ~ \hat{γ}_0 + \hat{γ}_1*κ to the sample {φ_h, h}^κ_{h=q}
          phismod <- log(1-phis); kappa <- length(InterKnots)
          gamma <- .lm.fit(cbind(1, q_F_X:kappa), phismod)$coef
          # Calculate \hat{φ}_κ based on the estimated coefficients
          phi_kappa <- 1 - exp(gamma[1])*exp(gamma[2]*kappa)
          # Store \hat{φ}_κ and the estimated coefficients \hat{γ}_0 and \hat{γ}_1
          phis_star <- c(phis_star, phi_kappa)
          oldintc <- c(oldintc, gamma[1]); oldslp <- c(oldslp, gamma[2])
          # Check if \hat{φ}_κ ≥ φ_{exit}
          if (phi_kappa >= phi_F_X) {
            #cat("Stopping iterations due to small improvement in SR\n\n")
            break
          }
        }
      }
    }



  }

  ########################
  # Final residuals plot #
  ########################
  if (plot && resids_plot && stop_plotting == max_iterations && stop_plotting == iter) {
    plot(X, resid_X, xlab = "x", ylab = "residuals", ylim = range(resid_X_1))
    legend("topleft",
           legend = c(expression(rho == F[N](x) - hat(F)[DDFS](x))),
           col = c("black"),
           pch = c(20, NA),
           bty = "n")
  }

  if (iter != max_iterations) {
    index <- iter-q_F_X
  } else {
    index <- iter
  }

  # For compatibility purposes (e.g. with GeDS::PPolyRep)
  newdata <- utils::getFromNamespace("read.formula", "GeDS")(formula = as.formula(F_X ~ f(X)), data = data.frame(X = X, F_X = F_X))
  F_X_hat_list[[index]]$terms <- newdata$terms

  out$f_X_hat <- f_X_hat_list[[index]]
  out$F_X_hat <- F_X_hat_list[[index]]
  out$RSS$mindist <- RSS; out$RSS$GeDS <- RSS_GeDS
  out$model <- index

  return(out)
}


