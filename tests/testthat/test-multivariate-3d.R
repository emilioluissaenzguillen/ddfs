test_that("three-dimensional fits define normalized joint distributions", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE),
              "Requires the multidimensional development version of GeDS")
  set.seed(903)
  XYZ <- cbind(x = runif(160), y = runif(160, -2, 3), z = runif(160, 10, 12))
  limits <- apply(XYZ, 2L, range)

  for (n in 2:4) {
    fit <- ddfs:::MultivariateDensityFitter(
      XYZ, n = n, min_iterations = 10L, max_iterations = 4L,
      tails_count_threshold = 0
    )
    expect_identical(fit$dimensions, colnames(XYZ))
    expect_equal(fit$ndim, 3L)
    expect_equal(fit$model, 4L)
    expect_equal(sum(lengths(fit$f_hat$knots) - 2L * n), 3L)
    expect_true(all(is.finite(fit$f_hat$coef)))
    expect_true(all(fit$f_hat$coef >= 0))

    # Evaluate the stored splines independently of the stored predictions.
    for (component in list(fit$f_hat, fit$F_hat)) {
      bases <- lapply(1:3, function(j) {
        splines::splineDesign(component$knots[[j]], XYZ[, j],
                              ord = component$order, outer.ok = TRUE)
      })
      values <- as.numeric(GeDS:::tensorProdND(bases) %*% component$coef)
      expect_equal(values, component$pred, tolerance = 1e-12)
    }

    grid <- as.matrix(expand.grid(lapply(1:3, function(j) {
      seq(limits[1L, j], limits[2L, j], length.out = 5L)
    })))
    bases <- lapply(1:3, function(j) {
      splines::splineDesign(fit$F_hat$knots[[j]], grid[, j],
                            ord = n + 1L, outer.ok = TRUE)
    })
    cdf <- as.numeric(GeDS:::tensorProdND(bases) %*% fit$F_hat$coef)
    expect_true(all(is.finite(cdf)))
    expect_true(all(cdf >= -1e-12 & cdf <= 1 + 1e-12))
    cdf_array <- array(cdf, c(5L, 5L, 5L))
    expect_equal(as.numeric(cdf_array[1L, , ]), rep(0, 25), tolerance = 1e-12)
    expect_equal(as.numeric(cdf_array[, 1L, ]), rep(0, 25), tolerance = 1e-12)
    expect_equal(as.numeric(cdf_array[, , 1L]), rep(0, 25), tolerance = 1e-12)
    expect_equal(cdf_array[5L, 5L, 5L], 1, tolerance = 1e-12)

    # A joint CDF must give non-negative probabilities to every grid cell,
    # not merely be non-decreasing in each coordinate separately.
    cells <- as.matrix(expand.grid(rep(list(1:4), 3L)))
    probabilities <- numeric(NROW(cells))
    for (mask in 0:7) {
      upper <- as.integer(intToBits(mask)[1:3])
      corners <- sweep(cells, 2L, upper, `+`)
      probabilities <- probabilities +
        (-1)^(3L - sum(upper)) * cdf_array[corners]
    }
    expect_true(all(probabilities >= -1e-12))
    expect_equal(sum(probabilities), 1, tolerance = 1e-12)

    # Independent integration oracle: two-point Gauss quadrature on every
    # knot interval is exact for the tested piecewise polynomial densities.
    queries <- rbind(
      limits[1L, ] + c(0.23, 0.61, 0.42) * (limits[2L, ] - limits[1L, ]),
      limits[1L, ] + c(0.78, 0.35, 0.86) * (limits[2L, ] - limits[1L, ]),
      limits[2L, ]
    )
    bases <- lapply(1:3, function(j) {
      splines::splineDesign(fit$F_hat$knots[[j]], queries[, j],
                            ord = n + 1L, outer.ok = TRUE)
    })
    cdf_at_queries <- as.numeric(GeDS:::tensorProdND(bases) %*% fit$F_hat$coef)
    for (i in seq_len(NROW(queries))) {
      axes <- lapply(1:3, function(j) {
        breaks <- sort(unique(c(fit$f_hat$knots[[j]][
          fit$f_hat$knots[[j]] < queries[i, j]], queries[i, j])))
        half_width <- diff(breaks) / 2
        midpoint <- head(breaks, -1L) + half_width
        cbind(nodes = c(midpoint - half_width / sqrt(3),
                        midpoint + half_width / sqrt(3)),
              weights = rep(half_width, 2L))
      })
      index <- as.matrix(expand.grid(lapply(axes, function(a) seq_len(NROW(a)))))
      points <- vapply(1:3, function(j) axes[[j]][index[, j], "nodes"],
                        numeric(NROW(index)))
      weights <- Reduce(`*`, lapply(1:3, function(j) {
        axes[[j]][index[, j], "weights"]
      }))
      bases <- lapply(1:3, function(j) {
        splines::splineDesign(fit$f_hat$knots[[j]], points[, j],
                              ord = n, outer.ok = TRUE)
      })
      density <- as.numeric(GeDS:::tensorProdND(bases) %*% fit$f_hat$coef)
      expect_true(all(density >= -1e-12))
      expect_equal(sum(weights * density), cdf_at_queries[i], tolerance = 1e-10)
      if (i == NROW(queries)) {
        expect_equal(sum(weights * density), 1, tolerance = 1e-10)
      }
    }
  }
})

test_that("three-dimensional tail constraints handle shared and repeated extrema", {
  skip_if_not(exists("tensorProdND", asNamespace("GeDS"), inherits = FALSE),
              "Requires the multidimensional development version of GeDS")
  set.seed(904)
  XYZ <- rbind(matrix(rnorm(600), ncol = 3), rep(-10, 3), rep(10, 3))
  colnames(XYZ) <- c("x", "y", "z")
  fit <- ddfs:::MultivariateDensityFitter(XYZ, n = 4L, max_iterations = 1L)
  repeated <- ddfs:::MultivariateDensityFitter(
    rbind(XYZ, XYZ[201:202, ]), n = 4L, max_iterations = 1L
  )
  expect_true(all(unlist(fit$args$tails)))
  index <- attr(fit$f_hat$coef, "basis.index")
  for (j in 1:3) {
    expect_true(all(fit$f_hat$coef[index[, j] %in% c(1L, 4L)] == 0))
  }
  expect_equal(fit$f_hat$pred[201:202], c(0, 0), tolerance = 1e-12)
  expect_equal(fit$F_hat$pred[201:202], c(0, 1), tolerance = 1e-12)
  expect_equal(as.numeric(fit$f_hat$coef), as.numeric(repeated$f_hat$coef),
               tolerance = 1e-10)
})
