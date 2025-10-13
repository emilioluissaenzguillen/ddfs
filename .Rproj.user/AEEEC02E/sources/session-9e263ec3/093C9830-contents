################################################################################
##################################### COEF #####################################
################################################################################
#' @title Coef method for ddfs objects
#' @name ddfs.GeDS
#' @description
#' Methods for the functions \code{\link[stats]{coef}} that allow to extract the
#' estimated coefficients of a fitted \code{\link{ddfs}}
#' object.
#' @param object the  \code{\link{ddfs}} object from which the
#' coefficients of the selected pdf/cdf should be extracted.
#' @param type Character string, either `"pdf"` or `"cdf"`, indicating which
#' coefficient vector to return; cdf coefficients are computed as weighted sums
#' of the pdf coefficients, where the weights depend on the knot locations.
#' @param ... Potentially further arguments (required by the definition of the
#' generic function). These will be ignored, but with a warning.
#'
#' @return A named vector containing the required coefficients of the fitted
#' univariate or bivariate ddfs model.
#'
#' @details
#' Simple method for the function \code{\link[stats]{coef}}.
#'
#' @seealso \code{\link[stats]{coef}} for the standard definition;
#' \code{\link{ddfs}} for examples.
#'
#' @export
#'
#' @aliases coef.ddfs
#' @rdname coef

coef.ddfs <- function(object, type = c("pdf", "cdf"), ...)
{

  # Handle additional arguments
  if(!missing(...)) warning("Only 'object', 'type' arguments will be considered")

  type <- match.arg(type)

  # map object$Type → dimension suffix
  suffix <- switch(object$Type,
                   "Univ - DDFS" = "_X_hat",
                   "Biv - DDFS"  = "_XY_hat",
                   stop("Unknown object$Type: ", object$Type))

  # prefix f vs F
  prefix <- if (type == "pdf") "f" else "F"

  # e.g. "f_X_hat" or "F_XY_hat"
  slot_name <- paste0(prefix, suffix)

  # extract and return the coef
  return(object[[slot_name]]$coef)
}

################################################################################
#################################### PREDICT ###################################
################################################################################
#' @title Density, Distribution, Quantile, and Random Sampling Methods for ddfs Objects
#'
#' @description
#' These functions provide access to the estimated model from a fitted
#' \code{ddfs} object. They allow evaluating the estimated probability density
#' function (\code{d.ddfs}), cumulative distribution function (\code{p.ddfs}),
#' quantile function (\code{q.ddfs}), and generate random samples (\code{r.ddfs}).
#'
#' @param object A \code{ddfs} object, typically returned by \code{\link{ddfs}}.
#' @param x,q Points at which to evaluate the density or distribution function.
#' @param p Probabilities at which to evaluate the quantile function.
#' @param N Number of random values to generate.
#'
#' @return
#' A numeric vector:
#' \itemize{
#'   \item \code{d.ddfs}: estimated density values at the specified points.
#'   \item \code{p.ddfs}: estimated cumulative distribution values at the
#'   specified quantiles.
#'   \item \code{q.ddfs}: estimated quantiles corresponding to the specified
#'   probabilities.
#'   \item \code{r.ddfs}: \code{N} random samples drawn from the estimated
#'   distribution.
#' }
#'
#' @examples
#' set.seed(123)
#'
#' # Simulate data
#' N <- 1000
#' sim <- sim.dist(N = N, ex = "GEV Type III")
#' X <- sim$X
#'
#' # Fit a DDFS model
#' GmodDens <- ddfs(X, n = 4, min.intknots = 3, phi_F = 0.3, q_F = 1, beta = 0.1)
#'
#' # Evaluate density, cdf, and quantiles on new data
#' X_new <- sim.dist(N = 10, ex = "GEV Type III")$X
#' d.ddfs(GmodDens, X_new)
#' p.ddfs(GmodDens, X_new)
#' q.ddfs(GmodDens, p = c(0.25, 0.5, 0.75))
#'
#' # Generate random samples
#' r.ddfs(GmodDens, N = 10)
#'
#' @seealso \code{\link{ddfs}}
#'
#' @rdname ddfs_distribution_methods
#' @export
d.ddfs <- function(object, x) {
  # Class check
  if (!inherits(object, "ddfs")) {
    stop("The 'object' must be of class 'ddfs'.")
  }
  predict.ddfs(object, newdata = x, type = "density")
}

#' @rdname ddfs_distribution_methods
#' @export
p.ddfs <- function(object, q) {
  # Class check
  if (!inherits(object, "ddfs")) {
    stop("The 'object' must be of class 'ddfs'.")
  }
  predict.ddfs(object, newdata = q, type = "distribution")
}

#' @rdname ddfs_distribution_methods
#' @export
q.ddfs <- function(object, p) {
  # Class check
  if (!inherits(object, "ddfs")) {
    stop("The 'object' must be of class 'ddfs'.")
  }

  if(object$Type != "Univ - DDFS") {
    stop("Quantile prediction is only available for univariate ddfs fits.")
  }
  predict.ddfs(object, newdata = p, type = "quantile")
}

#' @rdname ddfs_distribution_methods
#' @export
r.ddfs <- function(object, N) {
  # Class check
  if (!inherits(object, "ddfs")) {
    stop("The 'object' must be of class 'ddfs'.")
  }
  q.ddfs(object, runif(N))
}

#' @method predict ddfs
#' @keywords internal
#' @importFrom splines polySpline splineDesign
#' @noRd
predict.ddfs <- function(object, newdata, type = "density",...)
  {

  ####################################################
  ################ UNIVARIATE DENSITY ################
  ####################################################
  if (object$Type == "Univ - DDFS") {

    if (!missing(newdata)) {
      if (NCOL(newdata) != 1L) {
        stop("'newdata' must contain exactly only one variable for univariate DDFS fits.")
      }
      X_new <- as.numeric(newdata)
      extr <- range(X_new)

    }

    if (type == "density") {
      if (missing(newdata)) return(as.numeric(object$f_X_hat$pred))
      object <- object$f_X_hat
      theta <- object$coef
      n <- object$order
      # Design matrix
      basisMatrix <- splineDesign(knots = object$knots,
                                  derivs = rep(0,length(X_new)), x = X_new,
                                  ord = n, outer.ok = T)
      pred <- basisMatrix %*% theta

    } else if (type == "distribution") {
      if (missing(newdata)) return(as.numeric(object$F_X_hat$pred))
      object <- object$F_X_hat
      theta <- object$coef
      n <- object$order
      # The integral of a spline of degree n is a spline of degree n + 1.
      basisMatrix <- splineDesign(knots = object$knots,
                                  derivs = rep(0,length(X_new)), x = X_new,
                                  ord = n, outer.ok = T)
      pred <- basisMatrix %*% theta
      # Replace predictions for elements of X that are out of bounds
      pred[X_new <= min(object$knots)] <- 0
      pred[X_new >= max(object$knots)] <- 1

    } else if (type == "quantile") {
      if (missing(newdata)) X_new <- as.numeric(object$F_X_hat$pred)
      # Check
      eps <- 1e-8
      if (any(X_new < -eps | X_new > 1 + eps)) {
        stop("Error: newdata must be in the range [0,1] for quantile calculations.")
      }
      object <- object$F_X_hat
      # (i) Piecewise polynomial representation
      kn <- object$knots
      cf <- object$coef
      n <- object$order
      newlist <- list(knots = kn, coefficients = cf, order = n)
      class(newlist) <- c("nbSpline", "bSpline", "spline")
      xname <- attr(object$terms,"specials")$f-1
      xname <- attr(object$terms,"term.labels")[xname]
      xname <- substr(xname,3,(nchar(xname)-1))
      yname <- rownames(attr(object$terms,"factors"))[1]
      fortmp <- paste0(yname," ~ ", xname)
      attr(newlist,"formula") <- as.formula(fortmp)
      ppoly <- polySpline(newlist)

      # (ii) Invert polynomial
      pred <- PPolyInv(ppoly, X_new)

    }

    ####################################################
    ################ BIVARIATE DENSITY #################
    ####################################################
  } else if (object$Type == "Biv - DDFS") {

    if (NCOL(newdata) != 2L) {
      stop("'newdata' must contain two variables for bivariate DDFS fits.")
    }

    XY_new <- newdata
    Xextr <- range(XY_new[,1])
    Yextr <- range(XY_new[,2])

    if (type == "density") {
      if (missing(newdata)) return(as.numeric(object$f_XY_hat$pred))
      object <- object$f_XY_hat
      theta <- object$coef
      n <- object$order
      kntX <- object$knots$Xk
      kntY <- object$knots$Yk
      # Create spline basis matrix using specified knots, evaluation points and order
      basisMatrixX <- splineDesign(knots = kntX, x = XY_new[,1], ord = n,
                                   derivs = rep(0,length(XY_new[,1])), outer.ok = T)
      basisMatrixY <- splineDesign(knots = kntY, x = XY_new[,2], ord = n,
                                   derivs = rep(0,length(XY_new[,2])), outer.ok = T)

      basisMatrixbiv <- utils::getFromNamespace("tensorProd", "GeDS")(basisMatrixX, basisMatrixY)

      pred <- basisMatrixbiv %*% theta

      # Define bounds for X and Y
      X_min <- min(kntX)
      X_max <- max(kntX)
      Y_min <- min(kntY)
      Y_max <- max(kntY)
      # Identify out-of-bound indices for X and Y
      out_of_bounds_X_low <- XY_new[,1] <= X_min
      out_of_bounds_X_high <- XY_new[,1] >= X_max
      out_of_bounds_Y_low <- XY_new[,2] <= Y_min
      out_of_bounds_Y_high <- XY_new[,2] >= Y_max
      # Set predictions for out-of-bound values
      pred[out_of_bounds_X_low | out_of_bounds_Y_low] <- 0
      pred[out_of_bounds_X_high | out_of_bounds_Y_high] <- 0

    } else if (type == "distribution") {
      if (missing(newdata)) return(as.numeric(object$F_XY_hat$pred))
      object <- object$F_XY_hat
      theta <- object$coef
      n <- object$order
      kntX <- object$knots$Xk
      kntY <- object$knots$Yk

      # Replace out-of-bounds values with the nearest valid knot to avoid getting 0
      XY_new[,1] <- pmax(pmin(XY_new[,1], max(kntX)), min(kntX))
      XY_new[,2] <- pmax(pmin(XY_new[,2], max(kntY)), min(kntY))

      # Create spline basis matrix using specified knots, evaluation points and order
      basisMatrixX <- splineDesign(knots = kntX, x = XY_new[,1], ord = n,
                                   derivs = rep(0,length(XY_new[,1])), outer.ok = T)
      basisMatrixY <- splineDesign(knots = kntY, x = XY_new[,2], ord = n,
                                   derivs = rep(0,length(XY_new[,2])), outer.ok = T)

      basisMatrixbiv <- utils::getFromNamespace("tensorProd", "GeDS")(basisMatrixX, basisMatrixY)

      pred <- basisMatrixbiv %*% theta

      # # Define bounds for X and Y
      # X_min <- min(kntX)
      # X_max <- max(kntX)
      # Y_min <- min(kntY)
      # Y_max <- max(kntY)
      # # Identify out-of-bound indices for X and Y
      # out_of_bounds_X_low <- XY_new[,1] <= X_min
      # out_of_bounds_X_high <- XY_new[,1] >= X_max
      # out_of_bounds_Y_low <- XY_new[,2] <= Y_min
      # out_of_bounds_Y_high <- XY_new[,2] >= Y_max
      # # Set predictions for out-of-bound values
      # pred[out_of_bounds_X_low | out_of_bounds_Y_low] <- 0
      # pred[out_of_bounds_X_high & out_of_bounds_Y_high] <- 1


    }

  }


  return(as.numeric(pred))

}

################################################################################
##################################### KNOTS ####################################
################################################################################
#' @title Knots method for GeDS, GeDSboost, GeDSgam
#' @name knots.GeDS
#' @description
#' Method for the generic function \code{\link[stats]{knots}} that allows the
#' user to extract the vector of knots of a ddfs fit.
#' @param Fn the \code{\link{ddfs}} object from which the vector of knots
#' of the pdf or cdf spline model should be extracted.
#' @param type Character string, either `"pdf"` or `"cdf"`, indicating which knots
#' to return. Knots are identical except that the cdf, being one degree higher,
#' repeats the boundary knots once more. When `options = "internal"`, the two
#' sets of knots coincide exactly.
#' @param options a character string specifying whether "\code{all}" knots,
#' including the left-most and the right-most limits of the interval embedding
#' the observations (the default) or only the "\code{internal}" knots should be
#' extracted.
#' @param ... potentially further arguments (required for compatibility with the
#' definition of the generic function). Currently ignored, but with a warning.
#'
#' @return A vector in which each element represents a knot of the ddfs pdf/cdf
#' spline model.
#'
#' @details
#' This is a method for the function \code{\link[stats]{knots}} in the
#' \pkg{stats} package.
#'
#' @seealso \code{\link[stats]{knots}} for the definition of the generic function;
#' \code{\link{ddfs}} for examples.
#'
#' @export
#'
#' @aliases knots.ddfs
#' @rdname knots

knots.ddfs <- function(Fn, type = c("pdf", "cdf"), options = c("all","internal"), ...) {

  # Handle additional arguments
  if(!missing(...)) warning("Arguments other than 'Fn', 'type' and 'options' currenly igored. \n Please check if the input parameters have been correctly specified.")

  # validate args
  type    <- match.arg(type)
  options <- match.arg(options)

  # build slot name: prefix f/F and suffix _X_hat/_XY_hat
  prefix <- if (type == "pdf") "f" else "F"
  suffix <- switch(Fn$Type,
                   "Univ - DDFS" = "_X_hat",
                   "Biv - DDFS" = "_XY_hat",
                   stop("Unknown Fn$Type: ", Fn$Type))
  slot_name <- paste0(prefix, suffix)

  # extract knots
  kn <- Fn[[slot_name]]$knots

  # optionally reduce to internal knots using the same order
  if (options == "internal") {
    depth <- Fn[[slot_name]]$order
    kn    <- get_internal_knots(kn, depth = depth)
  }

  return(kn)

}

################################################################################
##################################### PLOT #####################################
################################################################################
#' @title Plot method for ddfs objects.
#' @name plot.ddfs
#' @description
#' Plot method for ddfs objects. Plots ddfs fits.
#' @param x a \code{\link{ddfs}} object from which the ddfs fit(s) should
#' be extracted.
#' @param f (optional) specifies the underlying function or generating process
#' to which the model was fit. This parameter is useful if the user wishes to
#' plot the specified function/process alongside the model fit and the data
#' @param main optional character string to be used as a title of the plot.
#' Alternatively, if a numeric vector is provided (e.g., containing true pdf/cdf
#' values), these values will be used to display a reference curve/dots on the plot.
#' @param legend.pos the position of the legend within the panel. See
#' \link[graphics]{legend} for details.
#' @param type character string specifying the type of plot required. Should be
#' set either to \code{"density"} or  \code{"distribution"}.
#' @param ... further arguments to be passed to the
#' \code{\link[graphics]{plot.default}} function.
#' @importFrom graphics plot points lines legend mtext abline
#' @importFrom plot3D persp3D text3D points3D segments3D
#' @export
#' @method plot ddfs

plot.ddfs <- function(x, f = NULL, main = NULL, legend.pos = NULL,
                      type = "density", ...)
  {

  # # Check if x is of class "ddfs"
  # if(!inherits(x, "ddfs")) {
  #   stop("The input 'x' must be of class 'ddfs'")
  # }

  # Other arguments passed to the function
  others <- list(...)

  if (x$Type == "Univ - DDFS") {
    X_new <- seq(min(x$f_X_hat$knots), max(x$f_X_hat$knots), diff(range(x$f_X_hat$knots))/1000)

    if(!is.null(f)) {
      if (is.function(f)) {
        f_values <- f(X_new)
      } else {
        stop("Error: 'f' (pdf/cdf) should be a function")
      }
    } else f_values <- NULL

    if (type == "density") {
      pred <- predict.ddfs(x, newdata = X_new, type = "density")
      ylab <- expression(f[X])
      if(!is.null(f)) {
        legend.text <- c(bquote(f(x)),
                         bquote(hat(f)[DDFS](x)))
        legend.col <- c("black", "red")
        pch = c(NA, NA); lwd = c(2, 2)
      } else {
        legend.text <- bquote(hat(f)[DDFS](x))
        legend.col <-  c("red")
        pch = c(NA); lwd = c(2)
      }
      if (is.null(legend.pos)) legend.pos <- "topright"

      knots_values <- round(x$f_X_hat$knots,2)
      intknt_text <- if (!is.null(get_internal_knots(knots_values, x$f_X_hat$order))) {
        paste0(", ", get_internal_knots(knots_values, x$f_X_hat$order))
      } else {
        ""
      }

      } else if (type == "distribution") {
        pred <- predict.ddfs(x, newdata = X_new, type = "distribution")
        ylab <- expression(F[X])
        if(!is.null(f)) {
          legend.text <- c(bquote(F[N](x)),
                           bquote(F(x)),
                           bquote(hat(F)[DDFS](x)))
          legend.col <- c("darkgrey", "black", "red")
          pch = c(20, NA, NA); lwd = c(NA, 2, 2)
          } else {
            legend.text <- c(bquote(F[N](x)), bquote(hat(F)[DDFS](x)))
            legend.col <- c("darkgrey", "red")
            pch = c(20, NA); lwd = c(NA, 2)
          }
        if (is.null(legend.pos)) legend.pos <- "bottomright"

        knots_values <- round(x$F_X_hat$knots,2)
        intknt_text <- if (!is.null(get_internal_knots(knots_values, x$F_X_hat$order))) {
          paste0(", ", get_internal_knots(knots_values, x$F_X_hat$order))
        } else {
          ""
        }
      }

    yylim <- range(c(pred, f_values))
    ylim <- if (!is.null(others$ylim)) others$ylim else yylim
    plot(X_new, pred, type = "n", col = "red", lwd = 2,
         ylab = ylab, xlab = "", main = main, ylim = ylim, cex.lab = 1.2)
    if (type == "distribution") points(x$Args$X, x$Args$ecdf, col = "darkgrey")
    lines(X_new, pred, col = "red", lwd = 2)
    if (!is.null(f)) lines(X_new, f_values, col = "black", lwd = 2)

    # Knots
    knt_text <- paste0(
      "a = ",
      knots_values[1],
      paste(intknt_text, collapse = ""),
      ", b = ",
      knots_values[length(knots_values)]
    )
    if (nchar(knt_text) > 40) {
      knt_text_split <- strsplit(knt_text, ",\\s*")[[1]]
      mid <- ceiling(length(knt_text_split) / 2)
      knt_text1 <- paste0(paste(knt_text_split[1:mid], collapse = ", "), ",")
      knt_text2 <- paste(knt_text_split[(mid + 1):length(knt_text_split)], collapse = ", ")

      knots_text1 <- if (type == "density") {
        bquote(italic(bold(t)[k * "," * n]) == italic(.(paste0("{", knt_text1))))
      } else if (type == "distribution") {
        bquote(italic(bold(t)[k * "," * n+1]) == italic(.(paste0("{", knt_text1))))
      }

      knots_text2 <- bquote(italic(.(paste0(knt_text2, "}"))))
    } else {
      knots_text1 <- if (type == "density") {
        bquote(italic(bold(t)[k * "," * n]) == italic(.(paste0("{", knt_text,"}"))))
      } else if (type == "distribution") {
        bquote(italic(bold(t)[k * "," * n+1]) == italic(.(paste0("{", knt_text,"}"))))
      }
      knots_text2 <- ""
    }

    mtext(knots_text1, side = 1, line = 3, cex = 0.9)
    mtext(knots_text2, side = 1, line = 4, cex = 0.9)

    for(knt in x$f_X_hat$knots) {
      abline(v = knt, col = "gray", lty = 2)
    }

    # Legend
    if (legend.pos != "none") {
      legend(legend.pos,
             legend = legend.text,
             col = legend.col,
             pch = pch, lwd = lwd,
             bty = "n")
    }


  } else if (x$Type == "Biv - DDFS") {

    N <- nrow(x$Args$XY)
    X <- x$Args$XY[,1]; Y <- x$Args$XY[,2]; ecdf <- x$Args$ecdf

    if (!is.null(f) && is.function(f)) {
      f_values <- f(x$Args$XY)
    } else {
      f_values <- f
    }

    X_new <- seq(from = min(X), to = max(X), length.out = 2*round(sqrt(N)))
    Y_new <- seq(from = min(Y), to = max(Y), length.out = 2*round(sqrt(N)))
    grid_matrix <- expand.grid(X = X_new, Y = Y_new)

    if (type == "density") {

      f_XY_hat_val <- predict.ddfs(x, newdata = grid_matrix, type = "density")
      f_XY_hat_mat <- matrix(f_XY_hat_val, nrow = 2*round(sqrt(N)))

      persp3D(x = X_new , y = Y_new, z = f_XY_hat_mat, phi = 25, theta = 50,
              xlim = range(X_new), ylim = range(Y_new), zlim = range(f_XY_hat_mat),
              ticktype = "detailed", expand = 0.7, colkey = FALSE, border = "black",
              xlab = "", ylab = "", zlab = " ")
      # Axis labels
      if (is.null(others$cex.lab)) cex.lab <- 1 else cex.lab <- others$cex.lab
      text3D(
        range(X_new)[1] - 0.2 * diff(range(X_new)),
        range(Y_new)[1] - 0.2 * diff(range(Y_new)),
        sum(range(f_XY_hat_mat))/2,
        labels = expression(f[X[1]]*","*phantom()[X[2]]),
        add = TRUE,
        cex = cex.lab
      )

      text3D(round(sum(range(X))/2), range(Y)[1] - 0.15*diff(range(Y)), 0,
             expression(X[1]), add = TRUE, cex = cex.lab, srt = -45)
      text3D(range(X)[2] + 0.15*diff(range(X)), -2*sum(range(Y))/2, 0,
             expression(X[2]), add = TRUE, cex = cex.lab, srt = 45)

      # Identify points where Z is greater/smaller than the model predicted value
      if(!is.null(f_values)) {
        tmp <- (x$f_XY_hat$pred - f_values > 0)
        points3D(x = X[!tmp], y = Y[!tmp], z = f_values[!tmp], col = "red", cex = 0.8, pch = 19, add = TRUE)
        points3D(x = X[tmp], y = Y[tmp], z = f_values[tmp], col = "blue", cex = 0.8, pch = 19, add = TRUE)
        legend("topleft", inset = c(0.2, 0.15),
               legend = c(expression(f >= hat(f)), expression(f < hat(f))),
               col = c("red", "blue"),
               pch = 19,
               bty = "n",
               cex = 1)
      }

    } else if (type == "distribution") {


      F_XY_hat_val <- predict.ddfs(x, newdata = grid_matrix, type = "distribution")
      F_XY_hat_mat <- matrix(F_XY_hat_val, nrow = 2*round(sqrt(N)))

      persp3D(x = X_new , y = Y_new, z = F_XY_hat_mat, phi = 25, theta = 50,
              xlim = range(X_new), ylim = range(Y_new), zlim = c(0, 1),
              ticktype = "detailed", expand = 0.7, colkey = FALSE, border = "black",
              xlab = "", ylab = "", zlab = " ")

      # Axis labels & title
      if (is.null(others$cex.lab)) cex.lab <- 1 else cex.lab <- others$cex.lab
      text3D(
        range(X_new)[1] - 0.2 * diff(range(X_new)),
        range(Y_new)[1] - 0.2 * diff(range(Y_new)),
        0.5,
        labels = expression(F[X[1]]*","*phantom()[X[2]]),
        add = TRUE,
        cex = cex.lab
      )
      text3D(round(sum(range(X))/2), range(Y)[1] - 0.15*diff(range(Y)), 0, expression(X[1]),
             add = TRUE, cex = cex.lab, srt = -45)
      text3D(range(X)[2] + 0.15*diff(range(X)), round(sum(range(Y))/2), 0, expression(X[2]),
             add = TRUE, cex = cex.lab, srt = 45)

      # Knots vectors
      # KnotsX
      knotsX_values <- round(x$F_XY_hat$knots$Xk,2)
      intkntX_text <- if (!is.null(get_internal_knots(knotsX_values,x$F_XY_hat$order))) {
        paste0(", ", get_internal_knots(knotsX_values,x$F_XY_hat$order))
      } else {
        ""
      }
      kntX_text <- paste0(
        "a = ",
        knotsX_values[1],
        paste(intkntX_text, collapse = ""),
        ", b = ",
        knotsX_values[length(knotsX_values)]
      )
      knotsX_text1 <- bquote(italic(bold(t)[1 * ";" * k[1] * "," * n[1]+1]) == italic(.(paste0("{", kntX_text,"}"))))
      knotsX_text2 <- ""

      if (is.null(others$cex.main)) cex.main <- 0.7 else cex.main <- others$cex.main

      mtext(knotsX_text1, side = 3, line = 2, adj = 0.5, cex = cex.main)
      mtext(knotsX_text2, side = 3, line = 1, adj = 0.5, cex = cex.main)
      # Knots Y
      knotsY_values <- round(x$F_XY_hat$knots$Yk,2)
      intkntY_text <- if (!is.null(get_internal_knots(knotsY_values,x$F_XY_hat$order))) {
        paste0(", ", get_internal_knots(knotsY_values,x$F_XY_hat$order))
      } else {
        ""
      }
      kntY_text <- paste0(
        "a = ",
        knotsY_values[1],
        paste(intkntY_text, collapse = ""),
        ", b = ",
        knotsY_values[length(knotsY_values)]
      )
      knotsY_text1 <- bquote(italic(bold(t)[2 * ";" * k[2] * "," * n[2]+1]) == italic(.(paste0("{", kntY_text,"}"))))
      knotsY_text2 <- ""
      mtext(knotsY_text1, side = 3, line = 0, adj = 0.5, cex = cex.main)
      mtext(knotsY_text2, side = 3, line = -1, adj = 0.5, cex = cex.main)

      # Knots segments
      Xextr <- range(X_new); Yextr <- range(Y_new)
      InterKnotsX <- x$F_XY_hat$knots$Xk
      InterKnotsY <- x$F_XY_hat$knots$Yk
      segments3D(InterKnotsX, rep(Yextr[1],length(InterKnotsX)), rep(0, length(InterKnotsX)),
                 InterKnotsX, rep(Yextr[1],length(InterKnotsX)), rep(0.5, length(InterKnotsX)),
                 lwd = 2, add = TRUE, col = "darkgrey")
      segments3D(rep(Xextr[2],length(InterKnotsY)), InterKnotsY, rep(0, length(InterKnotsY)),
                 rep(Xextr[2],length(InterKnotsY)), InterKnotsY, rep(0.5, length(InterKnotsY)),
                 lwd = 2, add = TRUE, col = "darkgrey")

      # Identify points where Z is greater/smaller than the model predicted value
      tmp <- (x$F_XY_hat$pred - ecdf > 0)
      points3D(x = X[!tmp], y = Y[!tmp], z = ecdf[!tmp], col = "red", cex = 0.8, pch = 19, add = TRUE)
      points3D(x = X[tmp], y = Y[tmp], z = ecdf[tmp], col = "blue", cex = 0.8, pch = 19, add = TRUE)
      legend("topleft", inset = c(0.2, 0.15),
             legend = c(expression(F[N] >= hat(F)), expression(F[N] < hat(F))),
             col = c("red", "blue"),
             pch = 19,
             bty = "n",
             cex = 1)



    }

  }


}

################################################################################
##################################### PRINT ####################################
################################################################################
#' @title Print method for ddfs objects
#' @name print.ddfs
#' @description
#' S3 method for the generic function \code{\link[base]{print}} that displays
#' a concise summary of a fitted \code{\link{ddfs}} model.
#' @param x a \code{ddfs}-class object, as returned by \code{\link{ddfs}}.
#' @param digits the number of significant digits to use when printing numeric summaries.
#' @param ... further arguments (ignored).
#' @return Invisibly returns the input object, augmenting it with a \code{Print} component
#'   containing the summary items.
#' @export
print.ddfs <- function(x,
                       digits = max(3L, getOption("digits") - 3L),
                       ...) {


  suffix <- switch(x$Type,
                   "Univ - DDFS" = "_X_hat",
                   "Biv - DDFS" = "_XY_hat",
                   stop("Unknown x$Type: ", x$Type))

  cat(paste0("\n", x$Type, ":\n"))

  ## Fucntion call
  cat("\nCall:\n", paste(deparse(x$extcall), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  ## Knots
  kn <- knots(x, options = "internal")
  cat(paste0("Number of internal knots: ", length(kn)))
  cat("\n")

  ## Spline orders (degree + 1)
  cat("\nSpline models orders (degree + 1):\n")
  cat("  pdf order: ", x[[paste0("f", suffix)]]$order, "\n", sep = "")
  cat("  cdf order: ", x[[paste0("F", suffix)]]$order, "\n\n", sep = "")

  ## RSS between empirical CDF and fitted CDF
  cat("RSS (empirical cdf vs. cdf model): ")
  rss_val <- x$RSS$mindist[[ x$model ]]
  cat(format(rss_val, digits = digits))
  cat("\n")


  invisible(x)
}


################################################################################
################################## ROUGHNESS ###################################
################################################################################
roughness_numeric <- function(x, f_hat, trim = 0.005) {


  if (trim != 0) {
    n <- length(x)
    tr <- round(trim*n)
    x <- x[(tr+1):(n-tr)]
    f_hat <- f_hat[(tr+1):(n-tr)]
  }
  # Input validation
  n <- length(x)
  if (n != length(f_hat)) stop("x and f_hat must be same length")
  if (n < 3) stop("Need at least 3 points to compute a second derivative.")

  # Check spacing
  dxs <- diff(x)
  h <- mean(dxs)  # Use average spacing
  if (max(abs(dxs - h)) > .Machine$double.eps^0.5 * max(abs(h))) {
    warning("Grid spacing not exactly equal: using average interval h = ", h)
  }

  # second-difference approximation to f'' at interior points
  f2 <- (f_hat[3:n] - 2*f_hat[2:(n-1)] + f_hat[1:(n-2)]) / h^2

  # approximate integral as sum f2^2 * h
  roughness <- sum(abs(f2)) * h
  return(roughness)
}


################################################################################
################################ RISK MEASURES #################################
################################################################################

## 1) VaR
#' Value at Risk (VaR) for \code{ddfs} objects
#'
#' Computes the \emph{Value at Risk} (quantile) at a given confidence level
#' for an object of class \code{ddfs}.
#'
#' @param object An object of class \code{ddfs}, typically returned from
#' \code{\link{ddfs}}.
#' @param alpha Numeric in (0,1). The confidence level for VaR, e.g.,
#' \code{alpha = 0.99} corresponds to the 99\% quantile.
#' @param ... Potentially further arguments (required by the definition of the
#' generic function). These will be ignored, but with a warning.
#'
#' @return A single \code{numeric} value giving the estimated quantile
#' (Value at Risk) at level \code{alpha}.
#' @examples
#' ## Danish Fire Loss Example
#' data("danish", package = "SMPracticals")
#' X <- sort(as.numeric(unlist(danish)))
#' hist(X)
#'
#' GmodDens <- ddfs(X, n = 4, min.intknots = 4, phi_F = 0.9, q_F = 2, beta = 0.6)
#'
#' alpha <- 0.99
#' VaR_empirical <- as.numeric(quantile(X, alpha))
#' TVaR_empirical <- mean(X[X >= VaR_empirical])
#' print(VaR_empirical); print(TVaR_empirical)
#'
#' VaR(GmodDens, alpha = alpha)
#' TVaR(GmodDens, alpha = alpha, tail = "right")
#'
#' @seealso \code{\link{TVaR}}
#' @export
VaR <- function(object, ...) {
  UseMethod("VaR")
}

#' @rdname VaR
#' @method VaR ddfs
#' @export
VaR.ddfs <- function(object, alpha, ...) {

  # Handle additional arguments
  if(!missing(...)) warning("Only 'object', 'type' arguments will be considered")

  # Class check
  if (!inherits(object, "ddfs")) {
    stop("The 'object' must be of class 'ddfs'.")
  }

  # alpha check
  if (missing(alpha)) {
    stop("Confidence level 'alpha' must be provided.")
  }
  return(as.numeric(predict.ddfs(object, newdata = alpha, type = "quantile")))
}


## 2) TVaR
#' Tail Value at Risk (TVaR) for \code{ddfs} objects
#'
#' Computes the \emph{Tail Value at Risk} (TVaR), at a given confidence level for
#' a fitted \code{ddfs} object.
#'
#' @param object An object of class \code{ddfs}, typically returned from
#' \code{\link{ddfs}}.
#' @param alpha Numeric value in (0,1) specifying the tail probability. For example,
#' if \code{alpha = 0.99} and \code{tail = "right"}, TVaR estimates the expected
#' loss given that the loss exceeds the 99th percentile.
#' @param tail A character string, either \code{"right"} or \code{"left"}, indicating
#' whether to compute the upper or lower tail expectation.
#' @param ... Potentially further arguments (required by the definition of the
#' generic function). These will be ignored, but with a warning.
#'
#' @return A single \code{numeric} value giving the estimated Tail Value at Risk
#' (i.e., conditional expectation in the tail above/below the \eqn{\alpha}-quantile).
#'
#' @details The right-tail TVaR at level \code{alpha} is defined as the expected loss
#' conditional on exceeding the \code{alpha}-quantile. The left-tail TVaR corresponds
#' to the expected value below the \code{alpha}-quantile.
#'
#' @seealso \code{\link{VaR}}
#' @export
TVaR <- function(object, ...) {
  UseMethod("TVaR")
}

#' @rdname TVaR
#' @method TVaR ddfs
#' @export
TVaR.ddfs <- function(object, alpha, tail = "right", ...) {

  # Handle additional arguments
  if(!missing(...)) warning("Only 'object', 'type' arguments will be considered")

  # Class check
  if (!inherits(object, "ddfs")) {
    stop("The 'object' must be of class 'ddfs'.")
  }
  # alpha check
  if (missing(alpha)) {
    stop("Confidence level 'alpha' must be provided.")
  }
  # alpha check
  if (missing(tail)) {
    stop("Distribution tail must be provided.")
  }
  theta <- object$f_X_hat$coef
  t <- object$f_X_hat$knots
  n <- object$f_X_hat$order
  VaR_alpha <- VaR.ddfs(object, alpha = alpha)
  p <- length(theta)
  result <- 0

  # Sum from i to p
  i_max <- p

  if (tail == "left") {
    for (i in 1:i_max) {
      # print(lhs_integral(i, n, t, VaR_alpha))
      # print(rhs_formula(i, n, t, VaR_alpha))
      result <- result + theta[i] * rhs_formula(i, n, t, VaR_alpha)
    }
    return(result/alpha)

    } else if (tail == "right") {
      for (i in 1:i_max) {
        # print(rhs_formula(i, n, t, round(t[i+n],8)))
        # print(lhs_integral(i, n, t, t[i+n]))
        #
        # print(rhs_formula(i, n, t, VaR_alpha))
        # print(lhs_integral(i, n, t, VaR_alpha))
        result <- result + theta[i] * (rhs_formula(i, n, t, t[i+n]) - rhs_formula(i, n, t, VaR_alpha))
      }
      return(result/(1-alpha))

    }
}


####################################
############ LEFT-TAIL #############
####################################
TVaR_num <- function(object, alpha) {
  VaR_gamma <- function(gamma) {
    VaR.ddfs(object, gamma)
  }
  res <- safe_integrate(VaR_gamma, lower = 0, upper = alpha)
  if (is.null(res)) return(NA_real_)
  return(res$value / alpha)
}

#' @importFrom splines splineDesign
TVaR_num_prime <- function(object, alpha) {

  theta <- object$f_X_hat$coef
  t <- object$f_X_hat$knots
  n <- object$f_X_hat$order
  VaR_alpha <- VaR.ddfs(object, alpha = alpha)

  f <- function(x) x * splineDesign(knots = t, x = x, ord = n, derivs = 0, outer.ok = TRUE) %*% theta


  res <- safe_integrate(f, lower = -Inf, upper = VaR_alpha)
  if (is.null(res)) return(NA_real_)
  return(res$value / alpha)
}

####################################
############ RIGHT-TAIL ############
####################################
TVaR_right_num <- function(object, alpha) {
  VaR_gamma <- function(gamma) {
    VaR.ddfs(object, gamma)
  }
  res <- safe_integrate(VaR_gamma, lower = alpha, upper = 1)
  if (is.null(res)) return(NA_real_)
  return(res$value / (1-alpha))
}

#' @importFrom splines splineDesign
TVaR_right_num_prime <- function(object, alpha) {

  theta <- object$f_X_hat$coef
  t <- object$f_X_hat$knots
  n <- object$f_X_hat$order
  VaR_alpha <- VaR.ddfs(object, alpha = alpha)

  f <- function(x) x * splineDesign(knots = t, x = x, ord = n, derivs = 0, outer.ok = TRUE) %*% theta


  res <- safe_integrate(f, lower = VaR_alpha, upper = Inf)
  if (is.null(res)) return(NA_real_)
  return(res$value / (1-alpha))
}

###################
##### HELPERS #####
###################
safe_integrate <- function(f, lower, upper, ..., max_subdiv = 10000L) {

  subdivs <- c(100L, 500L, 1000L, 2000L, 5000L, max_subdiv)

  for (s in subdivs) {
    res <- tryCatch(
      integrate(f, lower = lower, upper = upper, subdivisions = s, ...),
      error = function(e) NULL
    )
    if (!is.null(res)) return(res)
  }

  warning("Integration failed even with ", max_subdiv, " subdivisions.")
  return(NULL)
}

# 1) ∫_0^{x}  s * N_{i,n}(s)  ds  (numerical)
#' @importFrom splines splineDesign
lhs_integral <- function(i, n, knots, x) {
  left  <- knots[i]       # t_i
  right <- knots[i + n]   # t_{i+n}

  # 1) If x is at or below the left‐knot, integral is zero
  if (x <= left) return(0)

  lower <- left
  upper <- min(x, right)     # never integrate past support

  # 2) Otherwise integrate only over the slice where N_{i,n} is non-zero
  # # Less precise:
  # f <- function(s) {
  #   bs <- splineDesign(knots, s, ord = n, outer.ok = TRUE)
  #   s * bs[, i]
  # }
  # integrate(f, lower = lower, upper = upper)$value

  # More precise: break the interval at each interior knot
  inner <- knots[ knots > lower & knots < upper ]
  brks  <- c(lower, inner, upper)

  f <- function(s) {
    # vectorised evaluation
    bs <- splineDesign(knots, s, ord = n, outer.ok = TRUE)
    s  * bs[, i]
  }

  total <- 0
  for (k in seq_len(length(brks) - 1)) {
    total <- total +
      integrate(f,
                lower   = brks[k],
                upper   = brks[k + 1],
                rel.tol = 1e-8,
                abs.tol = 1e-8)$value
  }
  total

}

# 2) Right-hand side based on Bhatti–Bracken Proposition 1
#' @importFrom splines splineDesign
rhs_formula <- function(i, n, t, x) {

  i <- i + 1

  ik <- get_internal_knots(t, n)
  extr <- c(min(t), max(t))
  t <- sort(c(ik, rep(extr,n+1)))

  x <- min(x, max(t))

  basisMatrix <- splineDesign(knots = t, x = x, ord = n + 1, derivs = 0, outer.ok = TRUE)

  term1 <- x * sum(basisMatrix[, i:ncol(basisMatrix)])

  term2 <- 0
  for (j in i:ncol(basisMatrix)) {
    term2 <- term2 + I_0i(j, n + 1, t, x)
  }

  return( ((t[i+n] - t[i])/n) * (term1 - term2) )
}

#' @importFrom splines splineDesign
I_0i <- function(i, n, t, x){
  i <- i + 1
  ik <- get_internal_knots(t, n)
  extr <- c(min(t), max(t))
  t <- sort(c(ik, rep(extr,n+1)))

  x <- min(x, max(t))

  basisMatrix <- splineDesign(knots = t, x = x, ord = n + 1, derivs = 0, outer.ok = TRUE)

  return(((t[i+n] - t[i])/n)*sum(basisMatrix[, i:ncol(basisMatrix)]))
}

#' @importFrom splines splineDesign
I_0i_num <- function(i, n, t, x){

  f <- function(s) {
    splineDesign(knots = t, s, ord = n, outer.ok = TRUE)[,i]
  }

  return(integrate(f, lower = 0, upper = x)$value)

}

# # Alternatively
# quantile_function <- function(p, cdf_spline) {
#   uniroot(function(x) predict.ddfs(GmodDens, newdata = x, type = "distribution") - p, interval = range(GmodDens$Args$X))$root
# }
#
# # Example usage
# alpha = 0.95
# predict.ddfs(GmodDens, newdata = alpha, type = "quantile")
# quantile_function(alpha, GmodDens)
# quantile(X, alpha)


################################################################################
################################### SUMMARY ####################################
################################################################################
#' @title Summary method for ddfs
#' @name summary.ddfs
#' @description
#' Method for the generic function \code{\link[base]{summary}} that allows you to
#' print on screen the main information related to a fitted \code{\link{ddfs}}
#' model.
#' Similar to \code{\link{print.ddfs}} but with some extra detail.
#' @param object the \code{\link{ddfs}} object for which the main
#' information should be printed on screen.
#' @param ... potentially further arguments (required by the definition of the
#' generic function).
#' @seealso \code{\link{print.ddfs}}
#' @export
#' @aliases summary.ddfs

summary.ddfs <- function(object, ...)
{

  # 1) Core print
  print.ddfs(object, ...)

  # 2)
  cat("phi = ", object$Args$phi, "and q = ", object$Args$q, "(stopping rule parameters);\n")
  cat("beta = ", object$Args$beta)

  if (object$Args$beta == 0.5) {
    cat(", meaning that the within-cluster mean residual and the cluster range were considered equally important when placing the knots.\n")
  } else if (object$Args$beta > 0.5) {
    cat(", meaning that more weight was given to the within-cluster mean residual than to the cluster range when placing the knots.\n")
  } else if (object$Args$beta < 0.5) {
    cat(", meaning that more weight was given to the cluster range than to the within-cluster mean residual when placing the knots.\n")
  }


  invisible(object)
}
