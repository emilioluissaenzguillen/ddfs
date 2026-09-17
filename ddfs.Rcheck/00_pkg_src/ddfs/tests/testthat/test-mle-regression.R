# Frozen pre-optimization implementation, used only as a numerical reference.
original_constrMLE <- function(theta, indexes, basisMatrix, first_term, tol = 1e-6) {
  theta_prev <- theta
  repeat {
    theta_prev[indexes] <- theta[indexes]
    density <- basisMatrix %*% theta_prev
    density[density == 0] <- 1e-6
    theta <- first_term * theta_prev * colSums(basisMatrix / as.vector(density))
    if (max(abs(theta[indexes] - theta_prev[indexes])) < tol) break
  }
  theta
}

test_that("optimizer retains zero and tiny density handling", {
  basis <- rbind(c(0, 0), c(1, 0), c(.5, .5), c(0, 1))
  expect_equal(ddfs:::constrMLE_R(c(1, 1), 1:2, basis, c(1/3, 1/3)),
               original_constrMLE(c(1, 1), 1:2, basis, c(1/3, 1/3)), tolerance = 1e-12)
  expect_equal(ddfs:::constrMLE_R(1, 1L, matrix(1e-310), 1),
               original_constrMLE(1, 1L, matrix(1e-310), 1), tolerance = 1e-12)
})

test_that("full fits retain coefficients, knots and stopping decisions", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  set.seed(124)
  Sigma <- .5 ^ abs(outer(1:3, 1:3, "-"))
  gaussian <- MASS::mvrnorm(200, rep(0, 3), Sigma)
  mixture <- gaussian
  mixture[, 1] <- mixture[, 1] + rep(c(-2, 2), 100)
  for (X in list(gaussian[, 1], gaussian[, 1:2], gaussian, mixture, exp(gaussian))) {
    fit_once <- function() ddfs(X, n = 4, min.intknots = 10,
                                beta = .1, phi_F = .95, q_F = 2)
    reference <- local({
      testthat::local_mocked_bindings(constrMLE_R = original_constrMLE, .package = "ddfs")
      fit_once()
    })
    candidate <- fit_once()
    expect_equal(candidate$model, reference$model)
    expect_equal(length(candidate$RSS$mindist), length(reference$RSS$mindist))
    expect_equal(knots(candidate), knots(reference), tolerance = 1e-10)
    expect_equal(coef(candidate), coef(reference), tolerance = 1e-9)
    expect_equal(coef(candidate, fit = "cdf"), coef(reference, fit = "cdf"), tolerance = 1e-9)
    expect_equal(predict(candidate), predict(reference), tolerance = 1e-9)
    expect_equal(predict(candidate, fit = "cdf"), predict(reference, fit = "cdf"), tolerance = 1e-9)
    expect_equal(candidate$RSS, reference$RSS, tolerance = 1e-9)
  }
})
