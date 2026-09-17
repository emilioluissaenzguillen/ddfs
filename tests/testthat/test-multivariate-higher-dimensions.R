test_that("four- and five-dimensional public fits are normalized joint distributions", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE),
              "Requires the multidimensional development version of GeDS")
  for (dimension in 4:5) {
    set.seed(910 + dimension)
    observations <- matrix(runif(100 * dimension), ncol = dimension)
    observations <- sweep(sweep(observations, 2, seq_len(dimension), `*`),
                          2, -seq_len(dimension), `+`)
    colnames(observations) <- paste0("coordinate", seq_len(dimension))
    limits <- apply(observations, 2, range)
    for (order in 2:4) {
      # Three iterations, with stopping disabled until after the last iteration.
      fit <- ddfs(observations, n = order, min.intknots = 2, max.intknots = 2,
                  beta = .1, phi_F = .95, q_F = 1, tails_count_threshold = 0)
      expect_identical(fit$ndim, dimension)
      expect_identical(fit$dimensions, colnames(observations))
      expect_equal(fit$model, 3)
      expect_equal(sum(lengths(knots(fit, options = "internal"))), 2)
      expect_true(all(is.finite(coef(fit)) & coef(fit) >= 0))
      expect_equal(d.ddfs(fit, observations), predict(fit), tolerance = 1e-12)
      expect_equal(p.ddfs(fit, observations[, dimension:1]),
                   predict(fit, fit = "cdf"), tolerance = 1e-12)

      # Independent integration: two-point Gaussian quadrature is exact for
      # density orders 2:4 on each knot interval. Integrate each univariate
      # basis separately, then enumerate tensor coefficients in documented order.
      queries <- rbind(limits[1, ] + seq(.2, .8, length.out = dimension) *
                         (limits[2, ] - limits[1, ]), limits[2, ])
      colnames(queries) <- colnames(observations)
      coefficient_indices <- as.matrix(expand.grid(
        rev(lapply(lengths(fit$f_hat$knots) - order, seq_len))))
      coefficient_indices <- coefficient_indices[, dimension:1, drop = FALSE]
      for (row in 1:2) {
        integrated_bases <- lapply(seq_len(dimension), function(j) {
          breaks <- sort(unique(c(fit$f_hat$knots[[j]][
            fit$f_hat$knots[[j]] < queries[row, j]], queries[row, j])))
          widths <- diff(breaks) / 2
          midpoints <- head(breaks, -1) + widths
          nodes <- c(midpoints - widths / sqrt(3), midpoints + widths / sqrt(3))
          basis <- splines::splineDesign(fit$f_hat$knots[[j]], nodes, ord = order)
          colSums(basis * rep(widths, 2))
        })
        integrals <- Reduce(`*`, lapply(seq_len(dimension), function(j) {
          integrated_bases[[j]][coefficient_indices[, j]]
        }))
        mass <- sum(as.numeric(coef(fit)) * integrals)
        expect_equal(mass, p.ddfs(fit, queries[row, , drop = FALSE]), tolerance = 1e-10)
        if (row == 2) expect_equal(mass, 1, tolerance = 1e-10)
      }

      # Check all cells of a two-bin-per-axis partition by inclusion-exclusion.
      grid <- as.matrix(expand.grid(lapply(seq_len(dimension), function(j) {
        seq(limits[1, j], limits[2, j], length.out = 3)
      })))
      colnames(grid) <- colnames(observations)
      cdf <- p.ddfs(fit, grid)
      expect_true(all(is.finite(cdf) & cdf >= -1e-12 & cdf <= 1 + 1e-12))
      cdf_array <- array(cdf, rep(3, dimension))
      cells <- as.matrix(expand.grid(rep(list(1:2), dimension)))
      probabilities <- numeric(nrow(cells))
      for (mask in 0:(2^dimension - 1)) {
        upper <- as.integer(intToBits(mask)[seq_len(dimension)])
        corners <- sweep(cells, 2, upper, `+`)
        probabilities <- probabilities + (-1)^(dimension - sum(upper)) * cdf_array[corners]
      }
      expect_true(all(probabilities >= -1e-12))
      expect_equal(sum(probabilities), 1, tolerance = 1e-10)
      expect_equal(p.ddfs(fit, matrix(Inf, 1, dimension)), 1, tolerance = 1e-12)
      lower_faces <- matrix(Inf, dimension, dimension)
      diag(lower_faces) <- -Inf
      expect_equal(p.ddfs(fit, lower_faces), rep(0, dimension), tolerance = 1e-12)
      expect_equal(d.ddfs(fit, lower_faces), rep(0, dimension))
    }
  }
})

test_that("higher-dimensional coefficient limits cover the CDF and knot growth", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  for (dimension in 4:5) {
    set.seed(920 + dimension)
    observations <- matrix(runif(100 * dimension), ncol = dimension)
    # Initial quadratic PDF has 3^d coefficients, but its CDF needs 4^d.
    expect_error(ddfs:::MultivariateDensityFitter(
      observations, n = 3, max_iterations = 1, tails_count_threshold = 0,
      max.coef = 4^dimension - 1), "Tensor basis requires")
    fit <- ddfs:::MultivariateDensityFitter(
      observations, n = 3, max_iterations = 1, tails_count_threshold = 0,
      max.coef = 4^dimension)
    expect_length(fit$F_hat$coef, 4^dimension)
    expect_error(ddfs:::MultivariateDensityFitter(
      observations, n = 3, max_iterations = 2, tails_count_threshold = 0,
      max.coef = 4^dimension), "Tensor basis requires")
  }
})
