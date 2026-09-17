test_that("distribution methods reject objects of the wrong class", {
  expect_error(d.ddfs(list(), 0), "class 'ddfs'")
  expect_error(p.ddfs(list(), 0), "class 'ddfs'")
  expect_error(q.ddfs(list(), 0.5), "class 'ddfs'")
  expect_error(r.ddfs(list(), 1), "class 'ddfs'")
})

test_that("univariate density and distribution predictions respect boundaries", {
  fit <- mock_univariate_ddfs()

  expect_equal(d.ddfs(fit, c(0.25, 0.75)), c(1, 1))
  expect_equal(p.ddfs(fit, c(-1, 0, 0.5, 1, 2)), c(0, 0, 0.5, 1, 1))
  expect_equal(predict(fit, fit = "pdf"), fit$f_X_hat$pred)
  expect_equal(predict(fit, fit = "cdf"), fit$F_X_hat$pred)
})

test_that("method arguments are validated explicitly", {
  fit <- mock_univariate_ddfs()

  expect_error(predict(fit, 0.5, fit = "invalid"), "one of")
  expect_error(derive.ddfs(fit, 0.5, fit = "invalid"), "one of")
  expect_error(knots.ddfs(list()), "class 'ddfs'")
  expect_error(predict.ddfs(list(), 0.5), "class 'ddfs'")
})

test_that("bivariate stored predictions are available without newdata", {
  fit <- mock_bivariate_ddfs()

  expect_equal(predict(fit, fit = "pdf"), fit$f_XY_hat$pred)
  expect_equal(predict(fit, fit = "cdf"), fit$F_XY_hat$pred)
  expect_error(predict(fit, fit = "qf"), "only available for univariate")
})

test_that("coefficient and knot accessors select the requested fit", {
  fit <- mock_univariate_ddfs()

  expect_equal(coef(fit, fit = "pdf"), fit$f_X_hat$coef)
  expect_equal(coef(fit, fit = "cdf"), fit$F_X_hat$coef)
  expect_equal(knots(fit, fit = "pdf"), fit$f_X_hat$knots)
  expect_equal(knots(fit, fit = "cdf"), fit$F_X_hat$knots)
})
