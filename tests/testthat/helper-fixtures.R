mock_univariate_ddfs <- function() {
  knots_pdf <- c(rep(0, 2), rep(1, 2))
  knots_cdf <- c(rep(0, 3), rep(1, 3))

  structure(
    list(
      type = "Univ - DDFS",
      f_X_hat = list(pred = c(1, 1), knots = knots_pdf, coef = c(1, 1), order = 2L),
      F_X_hat = list(pred = c(0, 1), knots = knots_cdf, coef = c(0, 0.5, 1), order = 3L)
    ),
    class = "ddfs"
  )
}

mock_bivariate_ddfs <- function() {
  structure(
    list(
      type = "Biv - DDFS",
      f_XY_hat = list(pred = c(0.2, 0.3), coef = numeric(), order = 2L),
      F_XY_hat = list(pred = c(0.4, 0.8), coef = numeric(), order = 3L)
    ),
    class = "ddfs"
  )
}
