################################################################################
######################## Multivariate Density Fitter ###########################
################################################################################

# Fit a joint tensor-product DDFS model in two or more dimensions.
#'
#' @param data Numeric matrix or data frame with observations in rows and
#'   dimensions in columns.
#' @param n Spline order of the density estimate.
#' @param min_iterations Minimum number of knot-insertion iterations.
#' @param max_iterations Maximum number of knot-insertion iterations.
#' @param max.intknots Retained for compatibility. Each GeDS update inserts at
#'   most one new knot across all dimensions.
#' @param beta GeDS knot-placement parameter.
#' @param phi_F,q_F Parameters of the DDFS residual stopping rule.
#' @param tails_count_threshold Threshold used to identify decreasing tails on
#'   each coordinate.
#' @param plot,resids_plot,pdf,cdf,stop_plotting Plot controls. Direct plotting
#'   is currently available only for two-dimensional fits.
#' @param max.coef Maximum permitted tensor-basis size.
#' @return An internal DDFS fit with dimension-neutral `f_hat` and `F_hat`
#'   components.
#' @keywords internal
#' @noRd
#' @importFrom splines splineDesign
#' @importFrom utils getFromNamespace
MultivariateDensityFitter <- function(data, n = 3L, min_iterations = 1L,
                                      max_iterations = 50L, max.intknots = 1L,
                                      beta = 0.1, phi_F = 0.95, q_F = 2L,
                                      tails_count_threshold = 0.035,
                                      plot = FALSE, resids_plot = FALSE,
                                      pdf = NULL, cdf = NULL,
                                      stop_plotting = 0L,
                                      max.coef = 100000L)
{
  coordinates <- as.matrix(data)
  if (!is.numeric(coordinates)) {
    stop("Every dimension in 'data' must be numeric.", call. = FALSE)
  }
  storage.mode(coordinates) <- "double"
  ndim <- NCOL(coordinates)
  nobs <- NROW(coordinates)
  n <- as.integer(n)
  min_iterations <- as.integer(min_iterations)
  max_iterations <- as.integer(max_iterations)
  q_F <- as.integer(q_F)

  if (ndim < 2L) stop("'data' must contain at least two dimensions.", call. = FALSE)
  if (nobs < 1L || anyNA(coordinates) || any(!is.finite(coordinates))) {
    stop("'data' must contain finite, non-missing observations.", call. = FALSE)
  }
  if (any(vapply(seq_len(ndim), function(j) {
    length(unique(coordinates[, j])) < 2L
  }, logical(1)))) {
    stop("Every dimension in 'data' must contain at least two distinct values.",
         call. = FALSE)
  }
  if (is.null(colnames(coordinates))) colnames(coordinates) <- paste0("X", seq_len(ndim))
  if (anyDuplicated(colnames(coordinates))) {
    stop("Column names in 'data' must be unique.", call. = FALSE)
  }
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("'n' must be an integer greater than or equal to two.", call. = FALSE)
  }
  if (!is.numeric(max.coef) || length(max.coef) != 1L ||
      !is.finite(max.coef) || max.coef < 1) {
    stop("'max.coef' must be a positive number.", call. = FALSE)
  }

  dimension_names <- colnames(coordinates)
  coordinate_ranges <- lapply(seq_len(ndim), function(j) range(coordinates[, j]))
  names(coordinate_ranges) <- dimension_names
  empirical_cdf <- multivariate_ecdf(coordinates)
  residuals <- empirical_cdf
  internal_knots <- setNames(rep(list(NULL), ndim), dimension_names)
  density_fits <- cdf_fits <- list()
  RSS <- numeric()
  selected_index <- NA_integer_
  tail_constraints <- multivariate_tail_constraints(
    coordinates, threshold = tails_count_threshold
  )

  if (isTRUE(plot) && ndim > 2L) {
    warning("Direct plotting is currently available only for two-dimensional fits; fitting will continue without plots.")
    plot <- FALSE
  }

  for (iter in seq_len(max_iterations)) {
    if (iter > 1L) {
      if (isTRUE(plot) && isTRUE(resids_plot)) {
        # Retain the bivariate residual-fit display. This extra GeDS fit is
        # only for plotting; knot selection below always uses the ND engine.
        residual_plot <- withCallingHandlers(GeDS::NGeDS(
          residuals ~ f(X, Y),
          data = data.frame(X = coordinates[, 1L], Y = coordinates[, 2L]),
          beta = beta, max.intknots = iter - 1L,
          intknots_init = list(ikX = internal_knots[[1L]], ikY = internal_knots[[2L]]),
          higher_order = FALSE
        ), warning = function(w) {
          # This display fit is deliberately capped, not run to convergence.
          if (identical(conditionMessage(w), "Maximum number of iterations exceeded")) {
            invokeRestart("muffleWarning")
          }
        })
        plot(residual_plot, n = 2, xlab = "", ylab = "", zlab = "", main = "detail",
             legend.text = c(expression(rho >= hat(rho)), expression(rho < hat(rho))))
      }
      proposed_knots <- ddfs_insert_knot_nd(
        coordinates, residuals, internal_knots, coordinate_ranges,
        beta = beta, max.coef = max.coef
      )
      if (identical(proposed_knots, internal_knots)) {
        selected_index <- max(1L, iter - q_F)
        break
      }
      internal_knots <- proposed_knots
    }

    if (stop_plotting != 0L && stop_plotting + 1L == iter) {
      selected_index <- max(1L, iter - 1L)
      break
    }

    pdf_knots <- lapply(seq_len(ndim), function(j) {
      sort(c(internal_knots[[j]], rep(coordinate_ranges[[j]], n)))
    })
    names(pdf_knots) <- dimension_names
    basis_counts <- lengths(pdf_knots) - n
    # Check the larger CDF basis before allocating either dense design.
    basis_size <- prod(as.double(basis_counts + 1L))
    if (basis_size > max.coef) {
      stop("Tensor basis requires ", format(basis_size, scientific = FALSE),
           " coefficients, exceeding 'max.coef' = ",
           format(max.coef, scientific = FALSE), ".", call. = FALSE)
    }

    basis_matrices <- ddfs_basis_matrices(coordinates, pdf_knots, n)
    basis_matrix <- ddfs_tensor_product(basis_matrices)

    integral_widths <- Map(function(knots) diff(knots, lag = n) / n, pdf_knots)
    integral_weights <- ddfs_tensor_weights(integral_widths)
    basis_index <- ddfs_tensor_index(basis_counts, dimension_names)
    constrained <- ddfs_constrained_coefficients(
      basis_index, basis_counts, tail_constraints
    )
    active <- !constrained
    if (!any(active)) {
      stop("Tail constraints remove every tensor-product coefficient.", call. = FALSE)
    }

    theta <- numeric(length(integral_weights))
    theta[active] <- 1 / sum(integral_weights[active])
    # Constrained boundary faces have zero density, so observations on them
    # contribute nothing to the multiplicative MLE update. Count those rows
    # explicitly: subtracting the number of constrained faces (as in the
    # legacy fitter) is equivalent only when their extrema are distinct.
    effective_nobs <- sum(rowSums(basis_matrix[, active, drop = FALSE]) > 0)
    if (effective_nobs == 0L) {
      stop("No observations have positive support under the tail constraints.",
           call. = FALSE)
    }
    first_term <- 1 / (effective_nobs * integral_weights)
    theta <- constrMLE_R(theta, which(active), basis_matrix, first_term)
    theta[constrained] <- 0
    total_mass <- sum(theta * integral_weights)
    if (!is.finite(total_mass) || total_mass <= 0) {
      stop("Constrained density coefficients have invalid total mass.", call. = FALSE)
    }
    theta <- theta / total_mass
    density <- as.numeric(basis_matrix %*% theta)
    attr(theta, "basis.index") <- basis_index

    cdf_knots <- lapply(seq_len(ndim), function(j) {
      sort(c(internal_knots[[j]], rep(coordinate_ranges[[j]], n + 1L)))
    })
    names(cdf_knots) <- dimension_names
    cdf_basis_matrices <- ddfs_basis_matrices(coordinates, cdf_knots, n + 1L)
    cdf_basis <- ddfs_tensor_product(cdf_basis_matrices)
    theta_prime <- theta_prime_multivariate(theta, pdf_knots, rep(n, ndim))
    cdf_estimate <- as.numeric(cdf_basis %*% theta_prime)
    cdf_estimate <- pmin(pmax(cdf_estimate, 0), 1)
    residuals <- empirical_cdf - cdf_estimate
    RSS[iter] <- sum(residuals^2)

    density_fits[[iter]] <- list(
      pred = density, knots = pdf_knots, coef = theta, order = n,
      dimensions = dimension_names
    )
    cdf_fits[[iter]] <- list(
      pred = cdf_estimate, knots = cdf_knots, coef = theta_prime,
      order = n + 1L, dimensions = dimension_names
    )

    if (isTRUE(plot)) {
      temporary <- structure(list(
        f_XY_hat = ddfs_as_bivariate_component(density_fits[[iter]]),
        F_XY_hat = ddfs_as_bivariate_component(cdf_fits[[iter]]),
        type = "Biv - DDFS",
        args = list(XY = coordinates, ecdf = empirical_cdf)
      ), class = "ddfs")
      plot.ddfs(temporary, fit = "pdf", f = pdf)
      plot.ddfs(temporary, fit = "cdf", f = cdf)
    }

    if (iter < max_iterations && iter > q_F && iter > min_iterations &&
        RSS[iter] / RSS[iter - q_F] >= phi_F) {
      selected_index <- iter - q_F
      break
    }
  }

  if (!length(density_fits)) stop("No multivariate DDFS model was fitted.", call. = FALSE)
  if (isTRUE(plot) && isTRUE(resids_plot) &&
      stop_plotting == max_iterations && stop_plotting == iter) {
    plot3D::scatter3D(coordinates[, 1L], coordinates[, 2L], residuals,
                      xlab = "X", ylab = "Y", zlab = "residuals", colkey = FALSE)
  }
  if (is.na(selected_index)) selected_index <- length(density_fits)
  selected_index <- min(selected_index, length(density_fits))

  list(
    f_hat = density_fits[[selected_index]],
    F_hat = cdf_fits[[selected_index]],
    type = "Multiv - DDFS",
    dimensions = dimension_names,
    ndim = ndim,
    args = list(
      data = coordinates, ecdf = empirical_cdf, phi = phi_F, q = q_F,
      beta = beta, tails = tail_constraints, max.coef = max.coef
    ),
    RSS = list(mindist = RSS),
    model = selected_index
  )
}


# Preserve the established public bivariate component layout and plotting input.
ddfs_as_bivariate_component <- function(component)
{
  component$knots <- setNames(unname(component$knots), c("Xk", "Yk"))
  component$pred <- matrix(component$pred, ncol = 1L)
  component$coef <- as.numeric(component$coef)
  component$dimensions <- NULL
  component
}


# Request one additional knot from GeDS's dimension-independent Stage A.
ddfs_insert_knot_nd <- function(coordinates, residuals, internal_knots,
                                coordinate_ranges, beta, max.coef)
{
  make_grid <- utils::getFromNamespace("makeGridBoundsND", "GeDS")
  stage_a <- utils::getFromNamespace("stageALoopND", "GeDS")
  grid <- make_grid(coordinates, coordinate.ranges = coordinate_ranges)
  fit <- stage_a(
    coordinates = coordinates, response = residuals,
    upper.bounds = grid$upper.bounds, intknots = internal_knots,
    coordinate.ranges = coordinate_ranges,
    placement.ranges = coordinate_ranges, beta = beta,
    max.steps = 2L, stop.rule = "none", phi = 0.99, q = 1L,
    min.intknots = 0L, max.coef = max.coef
  )
  fit$selected.intknots
}


ddfs_basis_matrices <- function(coordinates, knots, order)
{
  matrices <- lapply(seq_len(NCOL(coordinates)), function(j) {
    splines::splineDesign(
      knots = knots[[j]], x = coordinates[, j], ord = order,
      derivs = rep(0L, NROW(coordinates)), outer.ok = TRUE
    )
  })
  names(matrices) <- colnames(coordinates)
  matrices
}


ddfs_tensor_product <- function(basis_matrices)
{
  utils::getFromNamespace("tensorProdND", "GeDS")(basis_matrices)
}


# Cartesian tensor index in GeDS column order: the last dimension varies fastest.
ddfs_tensor_index <- function(basis_counts, dimension_names = NULL)
{
  total <- prod(basis_counts)
  index <- vapply(seq_along(basis_counts), function(j) {
    block <- if (j == length(basis_counts)) {
      1
    } else {
      prod(basis_counts[(j + 1L):length(basis_counts)])
    }
    rep(rep(seq_len(basis_counts[j]), each = block), length.out = total)
  }, integer(total))
  colnames(index) <- if (is.null(dimension_names)) {
    paste0("Var", rev(seq_along(basis_counts)))
  } else {
    dimension_names
  }
  index
}


ddfs_tensor_weights <- function(weights)
{
  index <- ddfs_tensor_index(lengths(weights))
  apply(index, 1L, function(i) {
    prod(vapply(seq_along(weights), function(j) weights[[j]][i[j]], numeric(1)))
  })
}


ddfs_constrained_coefficients <- function(index, basis_counts, constraints)
{
  constrained <- rep(FALSE, NROW(index))
  for (j in seq_along(basis_counts)) {
    if (isTRUE(constraints[[j]]$lower)) constrained <- constrained | index[, j] == 1L
    if (isTRUE(constraints[[j]]$upper)) {
      constrained <- constrained | index[, j] == basis_counts[j]
    }
  }
  constrained
}


multivariate_tail_constraints <- function(data, threshold = 0.035)
{
  if (length(threshold) != 1L || is.na(threshold) || threshold < 0 || threshold > 1) {
    stop("'tails_count_threshold' must be a value in [0, 1].", call. = FALSE)
  }
  constraints <- lapply(seq_len(NCOL(data)), function(j) {
    values <- data[, j]
    limits <- range(values)
    width <- diff(limits)
    list(
      lower = mean(values <= limits[1L] + 0.05 * width) < threshold,
      upper = mean(values > limits[2L] - 0.05 * width) < threshold
    )
  })
  names(constraints) <- colnames(data)
  constraints
}
