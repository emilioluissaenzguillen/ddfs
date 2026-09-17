################################################################################
################################################################################
#################################### ddfs Class ################################
################################################################################
################################################################################
#' @title ddfs Class
#' @name ddfs-class
#' @description
#' A fitted ddfs object containing the results of a univariate or multivariate
#' ddfs model. This object stores both the estimated pdf and cdf models.
#'
#' Methods for functions \code{coef}, \code{knots}, \code{print}, \code{predict},
#' \code{plot}, as well as distribution-related methods such as \code{d.ddfs},
#' \code{p.ddfs}, \code{q.ddfs}, and \code{r.ddfs} are available.
#'
#' @section Components:
#' \describe{
#'   \item{f_hat, F_hat}{For fits above two dimensions, the pdf and cdf components,
#'   each containing predictions, a named knot list, coefficients, spline order,
#'   and dimension names. Bivariate components are named \code{f_XY_hat} and
#'   \code{F_XY_hat}; the univariate components are described below.}
#'   \item{dimensions, ndim}{For fits above two dimensions, the coordinate names
#'   and number of dimensions. See \code{\link{ddfs}} for supported methods.}
#'   \item{f_X_hat}{A list containing the estimated pdf model:
#' \itemize{
#'   \item pred a numeric matrix with the predicted density values at evaluation points.
#'   \item knots a numeric vector containing the locations of the knots in the spline fit.
#'   \item coef a numeric vector of estimated spline coefficients.
#'   \item order an integer representing the spline order used in the density estimation.
#'  }}
#'   \item{F_X_hat}{A list containing the estimated cdf model:
#' \itemize{
#'   \item pred a numeric vector containing the estimated cumulative probabilities.
#'   \item knots a numeric vector containing the locations of the knots in the cdf fit.
#'   \item coef a numeric vector of estimated spline coefficients.
#'   \item order an integer representing the spline order used in the CDF estimation.
#'  }}
#'   \item{type}{A character string specifying the type of estimation performed.
#' This can be \code{"Univ - DDFS"} for univariate ddfs estimation or
#'  \code{"Biv - DDFS"} for bivariate ddfs estimation, or
#'  \code{"Multiv - DDFS"} for higher-dimensional fits.}
#'
#'   \item{args}{A list of input parameters used in the estimation:
#' \itemize{
#'   \item phi A numeric value specifying a smoothing parameter.
#'   \item q A numeric value representing the penalty term in the estimation.
#'   \item beta A numeric value controlling the trade-off between smoothness and fit.
#'  }}
#'   \item{RSS}{A list containing residual sum of squares (RSS) measures:
#' \itemize{
#'   \item mindist A numeric vector of RSS values corresponding to the minimum distance method.
#'   \item GeDS A numeric vector of RSS values obtained using the GeDS approach.
#'  }}
#' }
#' @aliases ddfs-Class ddfs-class
#' @rdname ddfs-class
#'
#' @references
#' Dimitrova, D. S., Kaishev, V. K. and Saenz Guillen, E. (2025).
#' Distribution and density function estimation using variable-knot splines.
#' \emph{Manuscript submitted for publication.}
#' @docType class
NULL





