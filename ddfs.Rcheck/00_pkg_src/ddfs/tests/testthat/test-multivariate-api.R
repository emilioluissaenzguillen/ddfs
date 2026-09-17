test_that("public multidimensional fitting and methods agree with the fitter", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  set.seed(904)
  xyz <- cbind(x = runif(100), y = runif(100), z = runif(100))
  fit <- ddfs(xyz, n = 3, min.intknots = 1, max.intknots = 2,
              beta = .1, phi_F = .95, q_F = 1, tails_count_threshold = 0)
  reference <- ddfs:::MultivariateDensityFitter(
    xyz, n = 3, min_iterations = 2, max_iterations = 3,
    beta = .1, phi_F = .95, q_F = 1, tails_count_threshold = 0)
  expect_s3_class(fit, "ddfs")
  expect_identical(fit$type, "Multiv - DDFS")
  expect_equal(fit$f_hat, reference$f_hat)
  expect_equal(fit$F_hat, reference$F_hat)
  expect_equal(coef(fit), fit$f_hat$coef)
  expect_equal(coef(fit, fit = "cdf"), fit$F_hat$coef)
  expect_equal(knots(fit, options = "all"), fit$f_hat$knots)
  expect_equal(lengths(knots(fit, options = "internal")), lengths(fit$f_hat$knots) - 6L)
  expect_output(print(fit), "Dimensions: x, y, z")
  expect_output(summary(fit), "Multiv - DDFS")
  expect_equal(predict(fit), fit$f_hat$pred)
  expect_equal(predict(fit, fit = "cdf"), fit$F_hat$pred)
  expect_equal(d.ddfs(fit, xyz), fit$f_hat$pred)
  expect_equal(p.ddfs(fit, xyz), fit$F_hat$pred)
  expect_equal(d.ddfs(fit, as.data.frame(xyz[, 3:1])), fit$f_hat$pred)
  expect_equal(p.ddfs(fit, unname(xyz)), fit$F_hat$pred)
  expect_identical(d.ddfs(fit, xyz[FALSE, ]), numeric())
  expect_identical(p.ddfs(fit, xyz[FALSE, ]), numeric())

  corners <- rbind(c(Inf, Inf, Inf), c(-Inf, Inf, Inf), c(.5, Inf, Inf))
  expect_equal(p.ddfs(fit, corners)[1:2], c(1, 0), tolerance = 1e-12)
  expect_equal(d.ddfs(fit, corners), rep(0, 3))
  finite_corner <- matrix(c(.5, max(xyz[, 2]), max(xyz[, 3])), nrow = 1)
  expect_equal(p.ddfs(fit, corners)[3], p.ddfs(fit, finite_corner))
  expect_error(d.ddfs(fit, xyz[, 1:2]), "one column")
  expect_error(d.ddfs(fit, matrix("a", 1, 3)), "numeric")
  expect_error(d.ddfs(fit, matrix(NA_real_, 1, 3)), "missing")
  wrong <- xyz
  colnames(wrong) <- c("x", "y", "unknown")
  expect_error(d.ddfs(fit, wrong), "column names")
  colnames(wrong) <- c("x", "x", "z")
  expect_error(d.ddfs(fit, wrong), "column names")
  expect_error(q.ddfs(fit, .5), "univariate")
  expect_error(r.ddfs(fit, 2), "univariate")
  expect_error(derive.ddfs(fit, .5), "univariate")
  expect_error(plot(fit), "above two dimensions")
})

test_that("multidimensional public fitting makes current limits explicit", {
  xyz <- matrix(seq_len(60) / 60, ncol = 3)
  expect_error(ddfs(xyz, min.intknots = 1, beta = .1, phi_F = .95, q_F = 1,
                    stoptype = "SR"), "stoptype")
  expect_error(ddfs(xyz, min.intknots = 1, beta = .1, phi_F = .95, q_F = 1,
                    schoenberg = TRUE), "schoenberg")
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  expect_error(ddfs(xyz, n = 3, min.intknots = 1, max.intknots = 1,
                    beta = .1, phi_F = .95, q_F = 1, max.coef = 63),
               "Tensor basis requires")
  expect_error(ddfs(xyz, min.intknots = 1, beta = .1, phi_F = .95, q_F = 1,
                    max.coef = Inf), "max.coef")
})
