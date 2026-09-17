test_that("multivariate ECDF includes ties in every dimension", {
  x <- rbind(
    c(0, 0, 0),
    c(1, 0, 1),
    c(1, 1, 1)
  )

  expect_equal(ddfs:::multivariate_ecdf(x), c(1 / 3, 2 / 3, 1))
  expect_equal(
    ddfs:::multivariate_ecdf(x, rbind(c(-1, -1, -1), c(2, 2, 2))),
    c(0, 1)
  )
  expect_error(
    ddfs:::multivariate_ecdf(x, matrix(1, nrow = 1, ncol = 2)),
    "same number of dimensions"
  )
})

test_that("tensor indices follow GeDS coefficient ordering", {
  index <- ddfs:::ddfs_tensor_index(c(2, 3, 2), c("a", "b", "c"))

  expect_equal(NROW(index), 12L)
  expect_identical(index[1, ], c(a = 1L, b = 1L, c = 1L))
  expect_identical(index[2, ], c(a = 1L, b = 1L, c = 2L))
  expect_identical(index[NROW(index), ], c(a = 2L, b = 3L, c = 2L))
})

test_that("multivariate coefficient integration reproduces bivariate result", {
  knots_x <- c(0, 0, 0.5, 1, 1)
  knots_y <- c(-1, -1, 0, 1, 1)
  theta <- seq_len(9)

  expected <- ddfs:::theta_prime_bivariate_func(
    theta, knots_x, knots_y, n1 = 2, n2 = 2
  )
  actual <- ddfs:::theta_prime_multivariate(
    theta, list(X1 = knots_x, X2 = knots_y), orders = 2
  )

  expect_equal(as.numeric(actual), expected)
})

test_that("three-dimensional uniform density integrates to a valid CDF", {
  knots <- rep(list(c(0, 0, 1, 1)), 3)
  names(knots) <- paste0("X", 1:3)

  cdf_coef <- ddfs:::theta_prime_multivariate(rep(1, 8), knots, orders = 2)

  expect_length(cdf_coef, 27)
  expect_equal(cdf_coef[1], 0)
  expect_equal(tail(cdf_coef, 1), 1)
  expect_true(all(diff(sort(unique(cdf_coef))) >= 0))
})

test_that("tail constraints select complete tensor faces", {
  index <- ddfs:::ddfs_tensor_index(c(2, 3, 2))
  constraints <- list(
    list(lower = TRUE, upper = FALSE),
    list(lower = FALSE, upper = TRUE),
    list(lower = FALSE, upper = FALSE)
  )

  constrained <- ddfs:::ddfs_constrained_coefficients(
    index, c(2, 3, 2), constraints
  )

  expect_identical(constrained, index[, 1] == 1L | index[, 2] == 3L)
})
