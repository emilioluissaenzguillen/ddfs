################################################################################
################################################################################
##################################### DDFS #####################################
################################################################################
################################################################################
#' @title Density & Distribution Function variable-knot Spline estimation
#' @name ddfs
#' @description
#' The \code{ddfs()} function implements a non-parametric method for simultaneously
#' estimating the probability density function (pdf) and cumulative distribution
#' function (cdf) of a random variable using variable-knot spline models. The
#' function supports both univariate and bivariate density estimation.
#'
#' @param data One or two dimensional matrix/data frame, the raw data.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree
#' \eqn{+ 1}) of the ddfs density estimate. Note that the corresponding ddfs
#' cdf estimate will be of order \eqn{n+1}.
#' @param min.intknots Minimum number of knot-insertion iterations, \eqn{k - q},
#' that is, of internal knots in the final \code{ddfs} fits.
#' @param max.intknots Maximum number of knot-insertion iterations, \eqn{k - q},
#' corresponding to the number of internal knots in the final \code{ddfs} fits.
#' @param beta Numeric parameter in the interval \eqn{[0,1]} tuning the knot
#' insertion scheme applied during the fit of the residual vector
#' \eqn{\rho_i = F_N(X_i) - \hat{F}_{\mathrm{DDFS}}(X_i)}, \eqn{i = 1, \dots, N}.
#' This follows Stage A of the Geometrically Designed Splines (GeDS)
#' procedure (see Kaishev et al., 2016). The default is \code{0.1}.
#' @param phi_F numeric parameter in the interval \eqn{(0,1)}
#' specifying the threshold for the ddfs iterations stopping rule. Default
#' is equal to \code{0.3}.
#' @param q_F Integer parameter that fine-tunes the stopping rule in the \code{ddfs}
#' fitting procedure. The default is \code{1L}.
#' @param stoptype A character string indicating the type of ddfs stopping rule
#' to be used. It should be either one of \code{"RDMD"}, \code{"SRMD"}, \code{"RD"}
#' or \code{"SR"}. See details.
#' @param tail_decay Character string controlling whether tail-decay constraints
#' are imposed in the univariate fit. Use \code{"auto"} to choose the constrained
#' tail(s) by a histogram-based rule, \code{"none"} to disable the constraints,
#' \code{"left"} or \code{"right"} to constrain only one tail, or \code{"both"}
#' to constrain both tails.
#' @param tails_count_threshold Numeric threshold used only for bivariate fits.
#' For each margin, the algorithm computes the proportion of observations falling
#' in the lower and upper 5% of the observed range. If that proportion is below
#' this threshold, the corresponding tail is treated as decreasing and the
#' associated boundary constraint is imposed. Ignored for univariate fits.
#' @param plot Logical indicating whether to plot or not the pdf/cdf fits at
#' each iteration. Default is \code{FALSE}.
#' @param pdf An optional function representing the true probability density function
#' \eqn{f(X)}. If provided, it will be evaluated over a uniform grid and displayed
#' alongside the estimated density.
#' @param cdf An optional function representing the true cumulative distribution
#' function \eqn{F(X)}. If provided, it will be evaluated over a uniform grid and
#' displayed alongside the estimated distribution function.
#' @param resids_plot logical variable, set to \code{TRUE} if you'd also like the
#' residuals \eqn{\rho_i}, \eqn{i=1,...,N} to be plotted.
#' @param stop_plotting integer, specifies the knot-insertion iteration at which
#' to stop. This parameter is useful for controlling the plotting process,
#' particularly if you want to export the plots at intermediate steps.
#' @param schoenberg If apply averaging knot location at each iteration.
#'
#' @details For a detailed description of the methodology, refer to
#' Dimitrova et al. (2025) (forthcoming).
#' Different stopping rules are available for the \code{ddfs} iterative fitting procedure:
#' \itemize{
#'   \item \strong{RMD} (\emph{Ratio of Minimum Distances}): stopping rule based on the ratio
#'   of consecutive residual sum of squares (RSS) between the empirical distribution function
#'   \eqn{\hat{F}_N(x)} and the estimated distribution function \eqn{\hat{F}(x)}; see equation (6).
#'
#'   \item \strong{SRMD} (\emph{Smoothed Ratio of Minimum Distances}): smoothed version of RMD,
#'   as described in Section 3.1 of Dimitrova et al. (2023).
#'
#'   \item \strong{RD} (\emph{Ratio of Deviances}): stopping rule based on the ratio of consecutive
#'   deviances of the residual fit at each iteration, where the residuals are defined as
#'   \eqn{\rho_i = F_N(X_i) - F_{\mathrm{DDFS}}(X_i)},
#'   for \eqn{i = 1, \dots, N}.
#'
#'   \item \strong{SR} (\emph{Smoothed Ratio}): smoothed version of RD, also described in
#'   Section 3.1 of Dimitrova et al. (2023).
#' }
#'
#' @return #' A fitted ddfs object containing the results of a univariate or bivariate
#' ddfs model. This object stores both the estimated pdf and cdf models.
#'
#' Methods for functions \code{coef}, \code{knots}, \code{print}, \code{predict},
#' \code{plot}, as well as distribution-related methods such as \code{d.ddfs},
#' \code{p.ddfs}, \code{q.ddfs}, and \code{r.ddfs} are available.
#'
#' \describe{
#'   \item{f_X_hat}{A list containing the estimated pdf model:
#'     \itemize{
#'       \item pred a numeric matrix with the predicted density values at evaluation points.
#'       \item knots a numeric vector containing the locations of the knots in the spline fit.
#'       \item coef a numeric vector of estimated spline coefficients.
#'       \item order an integer representing the spline order used in the density estimation.
#'     }
#'     }
#'   \item{F_X_hat}{A list containing the estimated cdf model:
#'     \itemize{
#'       \item pred a numeric vector containing the estimated cumulative probabilities.
#'       \item knots a numeric vector containing the locations of the knots in the cdf fit.
#'       \item coef a numeric matrix of estimated spline coefficients.
#'       \item order an integer representing the spline order used in the CDF estimation.
#'       \item Type a character string specifying that this component represents a cdf.
#'     }
#'     }
#'   \item{Type}{A character string specifying the type of estimation performed.
#'   This can be \code{"Univ - DDFS"} for univariate ddfs estimation or
#'   \code{"Biv - DDFS"} for bivariate ddfs estimation.}
#'
#'   \item{Args}{A list of input parameters used in the estimation:
#'     \itemize{
#'       \item phi A numeric value specifying a smoothing parameter.
#'       \item q A numeric value representing the penalty term in the estimation.
#'       \item beta A numeric value controlling the trade-off between smoothness and fit.
#'     }
#'     }
#'   \item{RSS}{A list containing residual sum of squares (RSS) measures:
#'     \itemize{
#'       \item mindist A numeric vector of RSS values corresponding to the minimum distance method.
#'       \item GeDS A numeric vector of RSS values obtained using the GeDS approach.
#'     }
#'     }
#' }
#'
#' @examples
#' set.seed(123)
#' N <- 500
#' X <- rnorm(N)
#'
#' f_X <- dnorm(X)
#' F_X <- pnorm(X)
#'
#' ddfs_fit <- ddfs(X, n = 4, phi_F = 0.3, q_F = 1, plot = TRUE,
#'                  pdf = function(x) dnorm(x, mean = 0, sd = 1),
#'                  cdf = function(x) pnorm(x, mean = 0, sd = 1),
#'                  resids_plot = TRUE)
#'
#' print(ddfs_fit)
#' summary(ddfs_fit)
#' coef(ddfs_fit, fit = "pdf")
#' coef(ddfs_fit, fit = "cdf")
#' knots(ddfs_fit, options = "internal")
#'
#'\dontrun{
#' # Bivariate example
#' set.seed(123)
#' N <- 1000
#' rho <- 0.5
#'
#' sim <- sim.dist(N, ex = "BivGauss_0.5")
#' XY <- sim$X
#' X <- XY[,1]; Y <- XY[,2]
#' plot(sim$X)
#' f_XY <- sim$f_X
#'
#' ddfs_fit <- ddfs(data = XY, n = 4, min.intknots = 9, phi_F = 0.4, q_F = 2, beta = 0.1)
#'
#' print(ddfs_fit)
#' summary(ddfs_fit)
#' coef(ddfs_fit, fit = "pdf")
#' coef(ddfs_fit, fit = "cdf")
#' knots(ddfs_fit, options = "internal")
#' }
#'
#' @export
#'
#' @references
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S. and Verrall, R.J. (2016).
#' Geometrically designed, variable knot regression splines.
#' \emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
#' DOI: \doi{10.1007/s00180-015-0621-7}
#'
#' Dimitrova, D. S., Kaishev, V. K., Lattuada, A. and Verrall, R. J.  (2023).
#' Geometrically designed variable knot splines in generalized (non-)linear
#' models.
#' \emph{Applied Mathematics and Computation}, \strong{436}. \cr
#' DOI: \doi{10.1016/j.amc.2022.127493}
#'
#' Dimitrova, D. S., Kaishev, V. K. and Saenz Guillen, E. (2025).
#' Distribution and density function estimation using variable-knot splines.
#' \emph{Manuscript submitted for publication.}

ddfs <- function(data, n = 4L, min.intknots = NULL, max.intknots = 40L,
                 beta = NULL, phi_F = NULL, q_F = NULL, stoptype = "RDMD",
                 tail_decay = c("auto", "none", "left", "right", "both"),
                 tails_count_threshold = 0.05,
                 plot = FALSE, pdf = NULL, cdf = NULL, resids_plot = FALSE,
                 stop_plotting = 0L, schoenberg = FALSE) {

  extcall <- match.call()
  # Input checks for data structure
  if (missing(data) || is.null(data)) {
    stop("'data' must be provided.")
  }
  if (!(NCOL(data) %in% c(1L, 2L))) {
    stop("'data' must have exactly one column (univariate) or two columns (bivariate).")
  }

  # Automatic parameter selection
  # Check for missing or NULL parameters
  params_needed <- list(
    min.intknots = min.intknots,
    phi_F        = phi_F,
    q_F          = q_F,
    beta         = beta
  )

  # Identify which are missing
  missing <- names(params_needed)[sapply(params_needed, function(x) is.null(x) || length(x) == 0)]

  # If any are missing, call choose_params() once
  if (length(missing) > 0) {
    params_auto <- suppressMessages(suppressWarnings(choose_params(X = data, type = "bw_top_mean_Q3")))

    # Only fill in the missing ones
    list2env(params_auto[missing], envir = environment())

    message("Automatic parameter selection applied for: ", paste(missing, collapse = ", "))
  }


  min_iterations <- min.intknots + q_F
  max_iterations <- max.intknots + q_F

  # Parameter checks
  if (!is.numeric(n) || length(n) != 1 || n != as.integer(n) || n < 1)
    stop("'n' must be a positive integer.")
  if (!is.numeric(min.intknots) || length(min.intknots) != 1 || min.intknots != as.integer(min.intknots))
    stop("'min.intknots' must be a positive integer.")
  if (!is.numeric(max.intknots) || length(max.intknots) != 1 || max.intknots != as.integer(max.intknots) || max.intknots < min.intknots)
    stop("'max.intknots' must be an integer => min.intknots.")

  if (!is.numeric(beta) || length(beta) != 1 || beta < 0 || beta > 1)
    stop("'beta' must be a numeric scalar in the closed interval [0, 1].")
  if (!is.numeric(phi_F) || length(phi_F) != 1 || phi_F <= 0 || phi_F >= 1)
    stop("'phi_F' must be a numeric scalar in the open interval (0, 1).")
  if (!is.numeric(q_F) || length(q_F) != 1 || q_F <= 0 || q_F %% 1 != 0)
    stop("'q_F' must be a positive integer.")



  if (NCOL(data) == 1) {
    fit <- UnivariateDensityFitter(X = data, n = n, min_iterations = min_iterations,
                                   max_iterations = max_iterations, max.intknots = 1,
                                   beta = beta, phi_F_X = phi_F, q_F_X = q_F, stoptype = stoptype,
                                   tail_decay = tail_decay,
                                   plot = plot, resids_plot = resids_plot, pdf = pdf, cdf = cdf, stop_plotting = stop_plotting,
                                   schoenberg = schoenberg)
  } else if (NCOL(data) == 2) {
    fit <- BivariateDensityFitter(XY = data, n = n, min_iterations = min_iterations,
                                  max_iterations = max_iterations, max.intknots = 1,
                                  beta = beta, phi_F_XY = phi_F, q_F_XY = q_F,
                                  tails_count_threshold = tails_count_threshold,
                                  plot = plot, resids_plot = resids_plot, pdf = pdf, cdf = cdf, stop_plotting = stop_plotting)
  }

  fit$extcall <- extcall
  class(fit) <- "ddfs"

  return(fit)
}


