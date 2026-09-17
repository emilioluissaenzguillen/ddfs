# Direct regression checks against the unchanged bivariate implementation.
# These require the development GeDS source, not just its version number:
# released and development builds can both identify themselves as 0.3.5.
test_that("one-sided tails and stopping decisions match the legacy fitter", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE),
              "Requires the multidimensional development version of GeDS")
  set.seed(321)
  XY <- cbind(X = rexp(200), Y = rexp(200))
  constraints <- ddfs:::multivariate_tail_constraints(XY)
  expect_false(constraints$X$lower)
  expect_true(constraints$X$upper)
  expect_false(constraints$Y$lower)
  expect_true(constraints$Y$upper)

  for (q in 1:2) {
    legacy <- ddfs:::BivariateDensityFitter(
      XY, n = 4L, min_iterations = 3L, max_iterations = 15L,
      phi_F_XY = 0.5, q_F_XY = q
    )
    candidate <- ddfs:::MultivariateDensityFitter(
      XY, n = 4L, min_iterations = 3L, max_iterations = 15L,
      phi_F = 0.5, q_F = q
    )

    expect_lt(length(legacy$RSS$mindist), 15L)
    expect_equal(candidate$model, legacy$model)
    expect_equal(candidate$RSS$mindist, legacy$RSS$mindist, tolerance = 1e-9)
    expect_equal(unname(candidate$f_hat$knots),
                 unname(legacy$f_XY_hat$knots), tolerance = 1e-10)
    expect_equal(as.numeric(candidate$f_hat$coef),
                 as.numeric(legacy$f_XY_hat$coef), tolerance = 1e-9)
    expect_equal(as.numeric(candidate$F_hat$coef),
                 as.numeric(legacy$F_XY_hat$coef), tolerance = 1e-9)
    expect_equal(candidate$f_hat$pred, as.numeric(legacy$f_XY_hat$pred),
                 tolerance = 1e-9)
    expect_equal(candidate$F_hat$pred, as.numeric(legacy$F_XY_hat$pred),
                 tolerance = 1e-9)
  }
})

test_that("multivariate fitting reproduces unconstrained bivariate iterations", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE),
              "Requires the multidimensional development version of GeDS")
  set.seed(123)
  XY <- cbind(X = runif(100), Y = runif(100))

  for (n in 2:4) {
    for (iterations in 1:5) {
      legacy <- ddfs:::BivariateDensityFitter(
        XY, n = n, min_iterations = 10L, max_iterations = iterations,
        tails_count_threshold = 0
      )
      candidate <- ddfs:::MultivariateDensityFitter(
        XY, n = n, min_iterations = 10L, max_iterations = iterations,
        tails_count_threshold = 0
      )

      expect_equal(unname(candidate$f_hat$knots),
                   unname(legacy$f_XY_hat$knots), tolerance = 1e-10)
      expect_equal(unname(candidate$F_hat$knots),
                   unname(legacy$F_XY_hat$knots), tolerance = 1e-10)
      expect_equal(as.numeric(candidate$f_hat$coef),
                   as.numeric(legacy$f_XY_hat$coef), tolerance = 1e-9)
      expect_equal(as.numeric(candidate$F_hat$coef),
                   as.numeric(legacy$F_XY_hat$coef), tolerance = 1e-9)
      expect_equal(candidate$f_hat$pred, as.numeric(legacy$f_XY_hat$pred),
                   tolerance = 1e-9)
      expect_equal(candidate$F_hat$pred, as.numeric(legacy$F_XY_hat$pred),
                   tolerance = 1e-9)
      expect_equal(candidate$RSS$mindist, legacy$RSS$mindist, tolerance = 1e-9)
      expect_equal(candidate$model, legacy$model)
    }
  }
})

test_that("multivariate fitting reproduces bivariate constrained tails", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE),
              "Requires the multidimensional development version of GeDS")
  set.seed(123)
  XY <- cbind(X = rnorm(200), Y = rnorm(200))
  expect_true(all(unlist(ddfs:::multivariate_tail_constraints(XY))))

  for (n in 3:4) {
    for (iterations in 1:5) {
      legacy <- ddfs:::BivariateDensityFitter(
        XY, n = n, min_iterations = 10L, max_iterations = iterations
      )
      candidate <- ddfs:::MultivariateDensityFitter(
        XY, n = n, min_iterations = 10L, max_iterations = iterations
      )

      expect_equal(unname(candidate$f_hat$knots),
                   unname(legacy$f_XY_hat$knots), tolerance = 1e-10)
      expect_equal(as.numeric(candidate$f_hat$coef),
                   as.numeric(legacy$f_XY_hat$coef), tolerance = 1e-9)
      expect_equal(as.numeric(candidate$F_hat$coef),
                   as.numeric(legacy$F_XY_hat$coef), tolerance = 1e-9)
      expect_equal(candidate$f_hat$pred, as.numeric(legacy$f_XY_hat$pred),
                   tolerance = 1e-9)
      expect_equal(candidate$F_hat$pred, as.numeric(legacy$F_XY_hat$pred),
                   tolerance = 1e-9)
      expect_equal(candidate$RSS$mindist, legacy$RSS$mindist, tolerance = 1e-9)
      expect_equal(candidate$model, legacy$model)
    }
  }
})
