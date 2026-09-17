test_that("univariate boundary extension expands spline support", {
  set.seed(42)
  x <- rnorm(100)
  fit0 <- ddfs(x, min.intknots = 2, max.intknots = 3,
               beta = .1, phi_F = .5, q_F = 1,
               tail_decay = "both")
  fit1 <- ddfs(x, min.intknots = 2, max.intknots = 3,
               beta = .1, phi_F = .5, q_F = 1,
               tail_decay = "both", boundary_extension = 1)
  expect_equal(range(fit0$f_X_hat$knots), range(x))
  expect_lt(min(fit1$f_X_hat$knots), min(x))
  expect_gt(max(fit1$f_X_hat$knots), max(x))
  left <- ddfs(x, min.intknots = 2, max.intknots = 3,
               beta = .1, phi_F = .5, q_F = 1,
               tail_decay = "both", boundary_extension = c(1, 0))
  expect_lt(min(left$f_X_hat$knots), min(x))
  expect_equal(max(left$f_X_hat$knots), max(x))
  expect_error(ddfs(x, boundary_extension = -1), "non-negative")
})
