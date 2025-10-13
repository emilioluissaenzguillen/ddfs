################################################################################
################################################################################
#################################### ddfs Class ################################
################################################################################
################################################################################
#' @title ddfs Class
#' @name ddfs-class
#' @description
#' A fitted ddfs object containing the results of a univariate or bivariate
#' ddfs model. This object stores both the estimated pdf and cdf models.
#'
#' Methods for functions \code{coef}, \code{knots}, \code{print}, \code{predict},
#' \code{plot}, as well as distribution-related methods such as \code{d.ddfs},
#' \code{p.ddfs}, \code{q.ddfs}, and \code{r.ddfs} are available.
#'
#' @slot f_X_hat A list containing the estimated pdf model:
#' \itemize{
#'   \item pred a numeric matrix with the predicted density values at evaluation points.
#'   \item knots a numeric vector containing the locations of the knots in the spline fit.
#'   \item coef a numeric vector of estimated spline coefficients.
#'   \item order an integer representing the spline order used in the density estimation.
#'  }
#' @slot F_X_hat a list containing the estimated cdf model:
#' \itemize{
#'   \item pred a numeric vector containing the estimated cumulative probabilities.
#'   \item knots a numeric vector containing the locations of the knots in the cdf fit.
#'   \item coef a numeric matrix of estimated spline coefficients.
#'   \item order an integer representing the spline order used in the CDF estimation.
#'   \item Type a character string specifying that this component represents a cdf.
#'  }
#' @slot Type A character string specifying the type of estimation performed.
#' This can be \code{"Univ - DDFS"} for univariate ddfs estimation or
#'  \code{"Biv - DDFS"} for bivariate ddfs estimation.
#'
#' @slot Args A list of input parameters used in the estimation:
#' \itemize{
#'   \item phi A numeric value specifying a smoothing parameter.
#'   \item q A numeric value representing the penalty term in the estimation.
#'   \item beta A numeric value controlling the trade-off between smoothness and fit.
#'  }
#' @slot RSS A list containing residual sum of squares (RSS) measures:
#' \itemize{
#'   \item mindist A numeric vector of RSS values corresponding to the minimum distance method.
#'   \item GeDS A numeric vector of RSS values obtained using the GeDS approach.
#'  }
#' @aliases ddfs-Class ddfs-class
#' @rdname ddfs-class
#'
#' @references
#' Dimitrova, D. S., Kaishev, V. K. and Saenz Guillen, E. (2025).
#' Distribution and density function estimation using variable-knot splines.
#' \emph{Manuscript submitted for publication.}
#' @importFrom methods setClass

setClass(
  "ddfs",
  representation(
    f_X_hat = "list",
    F_X_hat = "list",
    Type = "character",
    Args = "list",
    RSS = "list"
  )
)





