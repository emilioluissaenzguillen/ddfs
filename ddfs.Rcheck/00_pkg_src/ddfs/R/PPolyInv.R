#####################
## Invert PPolyRep ##
#####################
#' @title Invert the piecewise polynomial representation of a spline object
#' @name PPolyInv
#' @description
#' Computes the inverse mapping of a piecewise polynomial spline object. Given a
#' strictly monotonic spline (produced by \code{\link[GeDS]{PPolyRep}} or similar),
#' the function returns the corresponding predictor values for a new set of
#' response values.
#'
#' @param ppoly A spline object of class \code{"npolySpline"}, \code{"polySpline"},
#' or \code{"spline"} that represents a piecewise polynomial form. The spline must be
#' strictly monotonic (either increasing or decreasing) to allow for inversion.
#' @param y_new A numeric vector of response values for which the corresponding predictor
#' values are sought.
#'
#' @return A numeric vector (or column matrix) of predictor values corresponding to the
#' input \code{y_new} values.
#'
#' @details
#' \code{PPolyInv} first verifies that the supplied \code{ppoly} object is
#' invertible by checking its strict monotonicity via the helper function
#' \code{is_invertible}. If the spline is not strictly monotonic, the function
#' stops with an error.
#'
#' The function extracts the knot locations and polynomial coefficients from
#' \code{ppoly} to build a data frame of polynomial segments. For each value in
#' \code{y_new}, it identifies the correct interval and uses the helper function
#' \code{solve_x} to solve the corresponding polynomial equation for the
#' predictor value.
#'
#' @examples
#' # Generate a data sample for the response variable
#' # Y and the single covariate X
#' library(GeDS)
#' set.seed(123)
#' N <- 1000
#' f_1 <- function(x) x^3
#' X <- sort(runif(N, min = -5, max = -3))
#' # Specify a model for the mean of Y to include only a component
#' # non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.2)
#'
#' Gmod <- NGeDS(Y ~ f(X), phi = 0.9)
#'
#' plot(Gmod)
#'
#' # Convert GeDS fit to a cubic piecewise polynomial representation
#' ppoly <- PPolyRep(Gmod, n = 4)
#' # Invert the spline using predicted values to recover predictor values
#' pred_new <- Gmod$cubic.fit$predicted
#' X_new <- PPolyInv(ppoly, pred_new)
#' # Compare recovered predictors to original values (differences should be near 0)
#' as.numeric(round(X_new - X, 4))
#'
#' @export
PPolyInv <- function(ppoly, y_new)
{
  # Check if ppoly is of one of the specified classes
  if (!inherits(ppoly, "npolySpline") && !inherits(ppoly, "polySpline") && !inherits(ppoly, "spline")) {
    stop("ppoly must be of class 'npolySpline', 'polySpline', or 'spline'")
  }
  # Check if ppoly is invertible
  is_invppoly <- is_invertible(ppoly)
  if(!is_invppoly$is_inv) stop("ppoly needs to be stricly monotonic to be invertible!")

  # Check if y_new is a numeric vector
  if (!is.numeric(y_new)) {
    stop("y_new must be a numeric vector")
  }

  # Ensure y_new is a column matrix
  y_new <- matrix(y_new, ncol = 1)
  # Track original order with id variable
  y_new <- cbind(y_new, id = 1:length(y_new))
  # Ensure y_new remains a matrix after ordering
  y_new <- y_new[order(y_new[, 1]), , drop = FALSE]

  # Knots & Coefficients #
  # 1) ppoly$knots: a vector of size k + 2 containing the complete set of
  # knots (internal knots plus the limits of the interval) of the GeDS fit.
  knots <- ppoly$knots; k <- length(knots) - 2
  # 2) Let us note that the first k (= number of internal knots) + 1 rows of the
  # matrix contain the n coefficients of the k + 1 consecutive pieces of the
  # piecewise polynomial representation.
  coefficients <- ppoly$coefficients[1:(k+1), , drop = FALSE]
  ylim <- ppoly$coefficients[,1]

  aux <- data.frame(
    start_x   = knots[1:(k+1)],
    end_x     = knots[2:(k+2)],
    start_y   = ylim[1:(k+1)],
    end_y     = ylim[2:(k+2)]
  )
  # Dynamically add polynomial terms based on the number of columns in `coefficients`
  poly_names <- c("constant", "linear", "quadratic", "cubic", "quartic",
                  "quintic", "sextic")
  for (i in 1:ncol(coefficients)) {
    term_name <- poly_names[i]
    aux[[term_name]] <- coefficients[, i]
  }


  x_new <- sapply(y_new[,1], solve_x, aux = aux,
                  slope = is_invppoly$slp, knots = knots, k = k)

  # x_new <- vector("numeric", length(y_new[,1]))
  #
  # for (i in 1:length(y_new[,1])) {
  #
  #   # First and last value
  #   if (abs(y_new[,1][i] - aux[1, "start_y"]) < 1e-5) {
  #     x_new[i] <- knots[1]
  #     next
  #   }
  #   if (abs(y_new[,1][i] - aux[nrow(aux), "end_y"]) < 1e-5) {
  #     x_new[i] <- knots[k+2]
  #     next
  #   }
  #
  #   interval_condition <- FALSE # Flag to track if the y interval is found
  #
  #   for (j in 1:NROW(aux)) {
  #
  #     if (is_invppoly$slp == "increasing") {
  #       interval_condition <- round(aux$start_y[j],5) <= round(y_new[,1][i],5) && y_new[,1][i] < aux$end_y[j]
  #     } else if (is_invppoly$slp == "decreasing") {
  #       interval_condition <- round(aux$start_y[j],5) >= round(y_new[,1][i],5) && y_new[,1][i] > aux$end_y[j]
  #     }
  #
  #     if (interval_condition) {
  #
  #       # p(x) = ax^2+bx+c = y ==> ax^2+bx+(c-y) = 0
  #       coef <- aux[,setdiff(names(aux), c("start_x", "end_x", "start_y", "end_y"))][j,]
  #       coef$constant <- coef$constant - y_new[,1][i]
  #
  #       ########################
  #       ######## Linear ########
  #       ########################
  #       if (length(coef) == 2) {
  #
  #         a <- coef$linear
  #         b <- coef$constant
  #         x_new[i]   <- -b/a + aux$start_x[j]
  #         break # exit the loop once the correct interval is found and processed
  #
  #       ########################
  #       ###### Quadratic #######
  #       ########################
  #       } else if (length(coef) == 3) {
  #
  #         a <- coef$quadratic
  #         b <- coef$linear
  #         c <- coef$constant
  #         discriminant <- if (is_invppoly$slp == "increasing") {
  #           sqrt(b^2 - 4*a*c)
  #           } else {
  #             -sqrt(b^2 - 4*a*c)
  #             }
  #         x_new[i]   <- (-b + discriminant) / (2*a) + aux$start_x[j]
  #         break # exit the loop once the correct interval is found and processed
  #
  #       ########################
  #       #### Cubic/Quartic #####
  #       ########################
  #       } else {
  #         roots <- polyroot(as.numeric(coef))
  #
  #         # Filter only real roots
  #         roots <- roots[abs(Im(roots)) < 1e-5]
  #         suppressWarnings({
  #           real_roots <- as.numeric(roots)
  #         })
  #         # If there is more than one root, choose the min non-neg
  #         if (length(real_roots) > 1) {
  #           if (all(real_roots < 0)) {
  #             real_roots <- real_roots[which.min(abs(real_roots))]
  #           } else {
  #             real_roots <- real_roots[real_roots >= 0]
  #             real_roots <- real_roots[which.min(real_roots)]
  #           }
  #         }
  #
  #         x_new[i]   <- real_roots + aux$start_x[j]
  #         break # exit the loop once the correct interval is found and processed
  #
  #       }
  #     }
  #   }
  # }

  # Recover original ordering based on "id" column
  y_new <- cbind(y_new, x_new)
  y_new <- y_new[order(y_new[, "id"]), , drop = FALSE]
  # Extract the "x_new" column as a column vector
  x_new <- y_new[, "x_new", drop = FALSE]

  return(x_new)
}

# Helper function to check if the PPoly is invertible; for continuous functions
# strict monotonicity is a necessary and sufficient condition.
is_invertible <- function(ps, grid_points = 1000) {
  # Define a grid over the domain
  xs <- seq(min(ps$knots), max(ps$knots), length.out = grid_points)
  # Predict the spline values on the grid
  ys <- predict(ps, xs)$y
  # Compute differences between consecutive values
  diffs <- diff(ys)
  # Check if the spline is strictly increasing or strictly decreasing
  is_increasing <- all(diffs > -1e-6)
  is_decreasing <- all(diffs < 1e-6)
  is_invertible <- is_increasing || is_decreasing

  # Return a list containing the invertibility and monotonicity status
  return(list(is_inv = is_invertible,
              slp = if (is_increasing) "increasing" else if (is_decreasing) "decreasing" else "not monotonic"))
}

# Helper to calculate polynomial roots
solve_x <- function(y_value, aux, slope, knots, k) {
  tolerance <- 1e-10

  # Return the domain endpoints exactly when the requested value is at
  # the boundary of the spline range.
  if (abs(y_value - aux$start_y[1]) <= tolerance) return(knots[1])
  if (abs(y_value - aux$end_y[nrow(aux)]) <= tolerance) return(knots[k + 2])

  # Locate a polynomial piece whose endpoint values contain y_value.
  lower_y <- pmin(aux$start_y, aux$end_y) - tolerance
  upper_y <- pmax(aux$start_y, aux$end_y) + tolerance
  candidates <- which(y_value >= lower_y & y_value <= upper_y)

  if (length(candidates) == 0L) {
    stop("The requested value lies outside the range of the spline.")
  }

  # At a shared knot either adjacent interval gives the same inverse.
  j <- candidates[1L]
  coefficient_names <- setdiff(
    names(aux),
    c("start_x", "end_x", "start_y", "end_y")
  )
  coefficients <- as.numeric(aux[j, coefficient_names, drop = TRUE])
  coefficients[1L] <- coefficients[1L] - y_value

  # polySpline coefficients use powers of the local coordinate
  # t = x - start_x. Solving only on the current interval prevents a
  # real root from another polynomial branch being selected.
  interval_width <- aux$end_x[j] - aux$start_x[j]
  polynomial <- function(t) {
    sum(coefficients * t^(seq_along(coefficients) - 1L))
  }

  left_value <- polynomial(0)
  right_value <- polynomial(interval_width)
  if (abs(left_value) <= tolerance) return(aux$start_x[j])
  if (abs(right_value) <= tolerance) return(aux$end_x[j])

  if (left_value * right_value > 0) {
    stop("The spline inverse could not be bracketed within its polynomial interval.")
  }

  local_root <- uniroot(
    polynomial,
    interval = c(0, interval_width),
    tol = .Machine$double.eps^0.5
  )$root

  aux$start_x[j] + local_root

}
