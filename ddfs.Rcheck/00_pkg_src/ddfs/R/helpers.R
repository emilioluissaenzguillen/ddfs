####################
## Bivariate ECDF ##
####################
bivariate_ecdf <- function(data, x_val, y_val) {
  n <- nrow(data)
  # Count the number of points (x_i, y_i) such that x_i <= x and y_i <= y
  count <- sum(data[, 1] <= x_val & data[, 2] <= y_val)
  # Compute the ECDF value
  ecdf_val <- count / n
  return(ecdf_val)
}

# Empirical joint CDF evaluated at arbitrary observations. Ties are included.
multivariate_ecdf <- function(data, newdata = data) {
  data <- as.matrix(data)
  newdata <- as.matrix(newdata)
  if (!is.numeric(data) || !is.numeric(newdata)) {
    stop("'data' and 'newdata' must be numeric.", call. = FALSE)
  }
  if (NCOL(data) != NCOL(newdata)) {
    stop("'data' and 'newdata' must have the same number of dimensions.",
         call. = FALSE)
  }
  vapply(seq_len(NROW(newdata)), function(i) {
    mean(rowSums(sweep(data, 2L, newdata[i, ], `<=`)) == NCOL(data))
  }, numeric(1))
}

################################################################################
########################## Integrate Bivariate Spline ##########################
################################################################################
#' @importFrom splines2 ibs
compute_bivariate_integral <- function(kntX, kntY, theta, n, X, Y) {
  stopifnot(length(X) == length(Y))

  # number of basis functions
  pX <- length(kntX) - n
  pY <- length(kntY) - n
  if (length(theta) != pX * pY) {
    stop(sprintf("theta must have length %d = %d*%d", pX*pY, pX, pY))
  }

  # reshape coefficients into matrix
  theta_mat <- matrix(theta, nrow = pX, ncol = pY, byrow = TRUE)

  # Boundary knots (needed explicitly in splines2)
  Boundary.knots.x <- c(kntX[n], kntX[length(kntX) - n + 1])
  Boundary.knots.y <- c(kntY[n], kntY[length(kntY) - n + 1])

  # Internal knots
  int_knots_x <- if ((length(kntX) - 2L * n) >= 1L) kntX[(n+1):(length(kntX)-n)] else NULL
  int_knots_y <- if ((length(kntY) - 2L * n) >= 1L) kntY[(n+1):(length(kntY)-n)] else NULL

  # Integrated basis matrices (rows = points, cols = basis fns)
  Ax <- ibs(x = X, degree = n - 1, knots = int_knots_x,
            Boundary.knots = Boundary.knots.x, intercept = TRUE)
  Ay <- ibs(x = Y, degree = n - 1, knots = int_knots_y,
            Boundary.knots = Boundary.knots.y, intercept = TRUE)

  Ax <- as.matrix(Ax)
  Ay <- as.matrix(Ay)

  # For each m: Ax[m,] %*% theta_mat %*% Ay[m,]^T
  vals <- rowSums((Ax %*% theta_mat) * Ay)
  as.numeric(vals)
}

# #' @importFrom reticulate r_to_py py_to_r
# compute_bivariate_integral <- function(kntX, kntY, theta, n, X, Y) {
#   # Degrees
#   kx <- as.integer(n - 1L)
#   ky <- as.integer(n - 1L)
#
#   # Ensure numpy arrays
#   tx <- r_to_py(as.numeric(kntX))
#   ty <- r_to_py(as.numeric(kntY))
#   c  <- r_to_py(as.numeric(theta))
#
#   # Lower bounds (as in LSQBivariateSpline.integral)
#   x0 <- as.numeric(kntX[n])
#   y0 <- as.numeric(kntY[n])
#
#   vals <- vapply(seq_along(X), function(i) {
#     xm <- X[i]; ym <- Y[i]
#     if (xm <= x0 || ym <= y0) return(0.0)
#     out <- scipy_fitpack$dblint(tx, ty, c, kx, ky, x0, xm, y0, ym)
#     # Explicit conversion to R
#     out_r <- py_to_r(out)
#     return(as.numeric(out_r))
#   }, 0.0)
#
#   return(vals)
# }

# compute_bivariate_integral_old <- function(X, Y, kntX, kntY, theta, n) {
#   # Initialize the bivariate spline
#   bivariate_spline <- scipy$LSQBivariateSpline(
#     x = X,
#     y = Y,
#     z = rep(1, length(X)),
#     tx = list(),
#     ty = list(),
#     kx = as.integer(n - 1),
#     ky = as.integer(n - 1)
#   )
#
#   # Set the spline knots and coefficients
#   bivariate_spline$tck[[1]] <- kntX
#   bivariate_spline$tck[[2]] <- kntY
#   bivariate_spline$tck[[3]] <- theta
#
#   # Define the integration limits
#   x_min <- kntX[1]
#   y_min <- kntY[1]
#
#   # Define the Python function using py_run_string()
#   reticulate::py_run_string("
# import numpy as np
# def compute_integrals(spline, x_min, y_min, X, Y):
#     X = np.array(X)
#     Y = np.array(Y)
#     return np.array([spline.integral(x_min, xm, y_min, ym) for xm, ym in zip(X, Y)])
# ")
#
#   # Call the Python function
#   F_XY_hat <- reticulate::py$compute_integrals(bivariate_spline, x_min, y_min, X, Y)
#
#   return(as.numeric(F_XY_hat))
# }


################################################################################
get_internal_knots <- getFromNamespace("get_internal_knots", "GeDS")


################################################################################
# Make sure cdf coefficients are strictly increasing
make_strict_inc <- function(x) {
  # 1) Cap to [0, 1] range
  x <- pmin(pmax(x, 0), 1)

  # 2) Enforce strict increase
  for (i in seq_len(length(x) - 1)) {
    if (x[i + 1] <= x[i]) {
      x[i + 1] <- x[i] + .Machine$double.eps
    }
  }

  x
}

################################################################################
#' @importFrom DescTools CCC KendallW
#' @importFrom diptest dip.test
#' @importFrom dccpp dcor
#' @importFrom entropy entropy discretize
#' @importFrom moments skewness kurtosis jarque.test
#' @importFrom nortest ad.test
#' @importFrom pcaPP cor.fk
#' @importFrom stats cov

features <- function(x) {

  if (NCOL(x) == 1) {
    # --- Univariate features ---
    skew     <- skewness(x)
    kurt     <- kurtosis(x)
    iqr2sd   <- IQR(x) / sd(x)
    rng2sd   <- diff(range(x)) / sd(x)
    cv       <- sd(x) / abs(mean(x))
    outliers <- sum(x < quantile(x, 0.25) - 1.5 * IQR(x) |
                      x > quantile(x, 0.75) + 1.5 * IQR(x))

    mad_val  <- mad(x)
    trimmed  <- mean(x, trim = 0.1)
    ent      <- entropy(discretize(x, numBins = 20))

    # P-value binning function
    bin_pval <- function(p) {
      cut(p,
          breaks = c(0, 0.01, 0.025, 0.05, 0.1, 1),
          labels = c("<=.01", ".01-.025", ".025-.05", ".05-.1", ">.1"),
          include.lowest = TRUE, right = FALSE)
    }

    # Categorized test p-values
    ad_pval_cat  <- bin_pval(ad.test(x)$p.value)
    # sp_pval_cat  <- bin_pval(shapiro.test(x)$p.value)
    dip_pval_cat <- bin_pval(dip.test(x)$p.value)
    jb_pval_cat  <- bin_pval(jarque.test(x)$p.value)

    # Final dataframe
    data.frame(
      skew, kurt, iqr2sd, rng2sd, cv, outliers,
      mad_val, trimmed, ent,
      ad_pval_cat, dip_pval_cat, jb_pval_cat
    )
  } else if (NCOL(x) == 2) {
    # --- Bivariate features ---
    X <- x[,1]
    Y <- x[,2]

    skewX <- skewness(X); skewY <- skewness(Y)
    kurtX <- kurtosis(X); kurtY <- kurtosis(Y)

    # Correlations
    pearson  <- cor(X, Y, method = "pearson")
    spearman <- cor(X, Y, method = "spearman")
    kendall  <- cor.fk(X, Y) # faster than cor(X, Y, method = "kendall")

    # Covariance matrix + determinant (generalized variance)
    covmat <- cov(cbind(X, Y))
    detcov <- det(covmat)


    # PCA-like features (orientation of scatter)
    eig <- eigen(covmat)
    major_axis <- eig$vectors[,1]
    eccen      <- sqrt(eig$values[1] / eig$values[2])  # eccentricity of ellipse

    # Distance correlation (nonlinear dependence)
    dcor_val <- dcor(X,Y) # faster than dcor(X, Y)

    # Concordance measures
    concord_Lin <- CCC(X,Y)$rho.c[[1]]

    # Outliers
    # Mahalanobis distance
    d2 <- mahalanobis(cbind(X, Y), colMeans(cbind(X, Y)), cov(cbind(X, Y)))
    # Threshold at chi-square quantile with df = 2
    cutoff <- qchisq(0.975, df = 2)
    out_biv <- sum(d2 > cutoff)


    data.frame(
      skewX, skewY, kurtX, kurtY,
      pearson, spearman, kendall,
      detcov, eccen, dcor_val,
      concord_Lin, out_biv
    )

  }

}

################################################################################

# Cache serialized selection models, not data-dependent predictions. Each loaded
# namespace gets its own cache; changed or removed files are checked on access.
.ddfs_parameter_models <- new.env(parent = emptyenv())

read_parameter_model <- function(path) {
  stamp <- file.info(path)[c("size", "mtime")]
  if (anyNA(stamp)) return(readRDS(path))
  cached <- .ddfs_parameter_models[[path]]
  if (is.null(cached) || !identical(cached$stamp, stamp)) {
    cached <- list(stamp = stamp, model = readRDS(path))
    .ddfs_parameter_models[[path]] <- cached
  }
  cached$model
}

choose_params <- function(X, type) {

  int1 <- 500
  int2 <- 5000

  if (NCOL(X) == 1) {
    N <- length(X)
    feat_X <- features(X)

    ## 1) min.intknots, min(k-q)

    if (startsWith(type, "bw")) {
      X_trim <- sort(X)
      # alpha
      if (N < int1) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_alpha_N=100_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_alpha_N=100_", type, ".rds"), package = "ddfs"))
      } else if (N >= int1 && N < int2) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_alpha_N=1000_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_alpha_N=1000_", type, ".rds"), package = "ddfs"))
      } else if (N >= int2) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_alpha_N=10000_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_alpha_N=10000_", type, ".rds"), package = "ddfs"))
      }
      alpha <- predict(Gmodboost, feat_X, n = 4)
      # The tuning experiment considered trimming proportions in [0, 0.10].
      # GeDSboost is an unconstrained regression model, so clip predictions to
      # that fitted parameter domain rather than allowing extrapolation.
      alpha <- clip_trimming_alpha(alpha)

      # bw
      if (N < int1) {
        # tree <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_bw_N=100_", type, ".rds"))
        tree <- read_parameter_model(system.file(paste0("choose_bw_N=100_", type, ".rds"), package = "ddfs"))
      } else if (N >= int1 && N < int2) {
        # tree <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_bw_N=1000_", type, ".rds"))
        tree <- read_parameter_model(system.file(paste0("choose_bw_N=1000_", type, ".rds"), package = "ddfs"))
      } else if (N >= int2) {
        # tree <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_bw_N=10000_", type, ".rds"))
        tree <- read_parameter_model(system.file(paste0("choose_bw_N=10000_", type, ".rds"), package = "ddfs"))
      }
      bw <- as.character(predict(tree, newdata = feat_X, type = "class"))

      if (feat_X$dip_pval_cat == ">.1") {
        if (N ==100) k_cap = 12 else k_cap = 30
        } else {
          k_cap = 35
        }

      min.intknots <- min_intknots_univ(X, alpha, method = bw, k_cap = k_cap)

    } else if (startsWith(type, "default")) {
      # min.intknots
      if (N < int1) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_min.intknots_N=100_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_min.intknots_N=100_", type, ".rds"), package = "ddfs"))
      } else if (N >= int1 && N < int2) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_min.intknots_N=1000_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_min.intknots_N=1000_", type, ".rds"), package = "ddfs"))
      } else if (N >= int2) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_min.intknots_N=10000_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_min.intknots_N=10000_", type, ".rds"), package = "ddfs"))
      }
      min.intknots <- round(predict(Gmodboost, feat_X, n = 4))

    }


    ## 2) phi_F and q_F
    # q_F
    if (N < int1) {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_q_F_N=100_", type, ".rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_q_F_N=100_", type, ".rds"), package = "ddfs"))
    } else if (N >= int1 && N < int2) {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_q_F_N=1000_", type, ".rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_q_F_N=1000_", type, ".rds"), package = "ddfs"))
    } else if (N >= int2) {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_q_F_N=10000_", type, ".rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_q_F_N=10000_", type, ".rds"), package = "ddfs"))
    }
    probs <- predict(Gmodboost, feat_X, n = 2)
    q_F <- ifelse(probs > 0.5, 2, 1)

    # phi_F
    if (q_F == 1) {
      if (N < int1) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_1_N=100_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_1_N=100_", type, ".rds"), package = "ddfs"))
      } else if (N >= int1 && N < int2) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_1_N=1000_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_1_N=1000_", type, ".rds"), package = "ddfs"))
      } else if (N >= int2) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_1_N=10000_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_1_N=10000_", type, ".rds"), package = "ddfs"))
      }

    } else if (q_F == 2) {
      if (N < int1) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_2_N=100_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_2_N=100_", type, ".rds"), package = "ddfs"))
      } else if (N >= int1 && N < int2) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_2_N=1000_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_2_N=1000_", type, ".rds"), package = "ddfs"))
      } else if (N >= int2) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_2_N=10000_", type, ".rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_2_N=10000_", type, ".rds"), package = "ddfs"))
      }
    }
    phi_F <- predict(Gmodboost, feat_X, n = 4)
    phi_F <- max(min(phi_F, 0.995), 0.1)  # ensure 0.1 <= phi_F <= 0.995

    ## 3) beta
    if (N < int1) {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_beta_N=100_", type, ".rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_beta_N=100_", type, ".rds"), package = "ddfs"))
    } else if (N >= int1 && N < int2) {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_beta_N=1000_", type, ".rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_beta_N=1000_", type, ".rds"), package = "ddfs"))
    } else if (N >= int2) {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_beta_N=10000_", type, ".rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_beta_N=10000_", type, ".rds"), package = "ddfs"))
    }

    beta <- predict(Gmodboost, feat_X, n = 4)
    beta <- min(1, max(0, beta)) # make sure it is in [0,1]

    return(list(min.intknots = min.intknots, beta = beta,
                phi_F = phi_F, q_F = q_F))

  } else if (NCOL(X) == 2) {

    type <- "bw_top_mean_Q3"

    N <- NROW(X)
    feat_X <- features(X)

    ## 1) min.intknots, min(k-q)
    # alpha
    alpha <- 0.025

    # bw
    if (N < int1) {
      # tree <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_bw_N=100_", type,"_biv.rds"))
      tree <- read_parameter_model(system.file(paste0("choose_bw_N=100_", type,"_biv.rds"), package = "ddfs"))
    } else {
      # tree <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_bw_N=1000_", type,"_biv.rds"))
      tree <- read_parameter_model(system.file(paste0("choose_bw_N=1000_", type,"_biv.rds"), package = "ddfs"))
    }
    bw <- as.character(predict(tree, newdata = feat_X, type = "class"))

    min.intknots <- min_intknots_biv(X, alpha, method = bw)

    ## 2) phi_F and q_F
    # q_F
    if (N < int1) {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_q_F_N=100_", type,"_biv.rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_q_F_N=100_", type,"_biv.rds"), package = "ddfs"))
    } else {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_q_F_N=1000_", type,"_biv.rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_q_F_N=1000_", type,"_biv.rds"), package = "ddfs"))
    }
    probs <- predict(Gmodboost, feat_X, n = 2)
    q_F <- ifelse(probs > 0.5, 2, 1)

    # phi_F
    if (q_F == 1) {
      if (N < int1) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_1_N=100_", type,"_biv.rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_1_N=100_", type,"_biv.rds"), package = "ddfs"))
      } else {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_1_N=1000_", type,"_biv.rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_1_N=1000_", type,"_biv.rds"), package = "ddfs"))
      }

    } else if (q_F == 2) {
      if (N < int1) {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_2_N=100_", type,"_biv.rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_2_N=100_", type,"_biv.rds"), package = "ddfs"))
      } else {
        # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_phi_F_2_N=1000_", type,"_biv.rds"))
        Gmodboost <- read_parameter_model(system.file(paste0("choose_phi_F_2_N=1000_", type,"_biv.rds"), package = "ddfs"))
      }
    }
    phi_F <- predict(Gmodboost, feat_X, n = 4)
    phi_F <- max(min(phi_F, 0.995), 0.1)  # ensure 0.1 <= phi_F <= 0.995

    ## 3) beta
    if (N < int1) {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_beta_N=100_", type,"_biv.rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_beta_N=100_", type,"_biv.rds"), package = "ddfs"))
    } else {
      # Gmodboost <- readRDS(paste0("/users/addj700/archive/ddfs/R/inst/choose_beta_N=1000_", type,"_biv.rds"))
      Gmodboost <- read_parameter_model(system.file(paste0("choose_beta_N=1000_", type,"_biv.rds"), package = "ddfs"))
    }

    beta <- predict(Gmodboost, feat_X, n = 4)
    beta <- min(1, max(0, beta)) # make sure it is in [0,1]

    return(list(min.intknots = min.intknots, beta = beta,
                phi_F = phi_F, q_F = q_F))

  }
}

################################################################################
#' @importFrom stats quantile
# robust range: 2.5%-97.5% trimmed range per axis
trim_tails <- function(v, p) {
  qs <- quantile(v, c(p, 1 - p), names = FALSE)
  v[v >= qs[1] & v <= qs[2]]
}

#' @importFrom stats bw.nrd0 bw.nrd bw.ucv bw.bcv bw.SJ
#' @importFrom moments kurtosis
min_intknots_univ <- function(X, alpha,
                              method = c("bw.nrd0", "bw.nrd", "bw.ucv", "bw.bcv", "bw.SJ"),
                              k_cap = 35) {

  stopifnot(is.vector(X))
  method <- match.arg(as.character(method),
                      c("bw.nrd0", "bw.nrd", "bw.ucv", "bw.bcv", "bw.SJ"))

  # Bandwidth
  bw <- do.call(method, list(X))

  # Robust ranges
  X_trimmed <- trim_tails(X, alpha)
  r <- diff(range(X_trimmed))

  # Internal knot counts
  k <- (r / bw) - 1

  # Kurtosis-based deflation
  c <- kurtosis(X_trimmed)

  k_raw <- k * (1/c)

  k_min <- round(min(k_raw, k_cap))
  k_min <- max(1L, k_min)   # at least 1 internal knots
  as.integer(k_min)

}

clip_trimming_alpha <- function(alpha, upper = 0.1) {
  min(max(as.numeric(alpha), 0), upper)
}


#' @importFrom ks Hns.diag Hpi.diag Hbcv.diag Hlscv.diag Hscv.diag
#' @importFrom moments kurtosis
min_intknots_biv <- function(X, alpha,
                             method = c("Hns.diag", "Hpi.diag", "Hbcv.diag", "Hlscv.diag", "Hscv.diag")) {

  stopifnot(is.matrix(X), ncol(X) == 2L)
  method <- match.arg(as.character(method),
                      c("Hns.diag","Hpi.diag","Hbcv.diag","Hlscv.diag","Hscv.diag"))

  H <- do.call(method, list(X))

  # Bandwidths along axes
  hX <- sqrt(H[1,1])
  hY <- sqrt(H[2,2])

  # Robust ranges
  X_trimmed <- trim_tails(X[,1], alpha); Y_trimmed <- trim_tails(X[,2], alpha)
  rX <- diff(range(X_trimmed)); rY <- diff(range(Y_trimmed))

  # Internal knot counts
  kX <- (rX / hX) - 1
  kY <- (rY / hY) - 1

  # Kurtosis-based deflation
  cx <- kurtosis(X_trimmed)
  cy <- kurtosis(Y_trimmed)

  k_raw <- kX * (1/cx) + kY * (1/cy)
  k_cap <- 35

  k_min <- round(min(k_raw, k_cap))
  k_min <- max(2L, k_min)   # at least 2 internal knots
  as.integer(k_min)

}
