test_that("higher-dimensional defaults fill missing parameters without automatic tuning", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  testthat::local_mocked_bindings(
    choose_params = function(...) stop("Automatic tuning must not be called"),
    MultivariateDensityFitter = function(...) list(received = list(...)),
    .package = "ddfs")
  for (dimension in 3:6) {
    X <- matrix(seq_len(10 * dimension), ncol = dimension)
    args <- ddfs(X)$received
    expect_equal(args$n, 4L)
    expect_equal(args$min_iterations, 12L)
    expect_equal(args$max_iterations, 42L)
    expect_equal(args$beta, .1)
    expect_equal(args$phi_F, .95)
    expect_equal(args$q_F, 2L)
  }
  X <- matrix(seq_len(30), ncol = 3)
  expect_equal(ddfs(X, min.intknots = NULL, beta = NULL,
                    phi_F = NULL, q_F = NULL)$received, ddfs(X)$received)
  args <- ddfs(X, n = 3, min.intknots = 2, beta = .3, q_F = 1)$received
  expect_equal(args$n, 3)
  expect_equal(args$min_iterations, 3)
  expect_equal(args$beta, .3)
  expect_equal(args$q_F, 1)
  expect_equal(args$phi_F, .95)
  expect_equal(ddfs(X, phi_F = .7)$received$phi_F, .7)
})

test_that("default and explicit higher-dimensional fits agree", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  set.seed(941)
  X <- matrix(rnorm(300), ncol = 3)
  automatic <- ddfs(X, max.intknots = 10)
  explicit <- ddfs(X, n = 4, min.intknots = 10, beta = .1,
                   phi_F = .95, q_F = 2, max.intknots = 10)
  automatic$extcall <- explicit$extcall <- NULL
  expect_equal(automatic, explicit)
})

test_that("univariate automatic selection still fills only missing values", {
  testthat::local_mocked_bindings(
    choose_params = function(...) list(min.intknots = 2, beta = .1, phi_F = .4, q_F = 1),
    UnivariateDensityFitter = function(...) list(received = list(...)),
    .package = "ddfs")
  expect_message(fit <- ddfs(1:20, beta = .3), "Automatic parameter selection")
  expect_equal(fit$received$min_iterations, 3)
  expect_equal(fit$received$beta, .3)
  expect_equal(fit$received$phi_F_X, .4)
  expect_equal(fit$received$q_F_X, 1)
})
