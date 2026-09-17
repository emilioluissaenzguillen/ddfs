test_that("public bivariate fits preserve the legacy layout and numerical results", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  set.seed(123)
  samples <- list(cbind(X = runif(100), Y = runif(100)),
                  cbind(X = rnorm(200), Y = rnorm(200)))
  for (sample in seq_along(samples)) {
    XY <- samples[[sample]]
    threshold <- if (sample == 1) 0 else .035
    for (order in (if (sample == 1) 2:4 else 3:4)) {
      for (q in 1:2) {
        legacy <- ddfs:::BivariateDensityFitter(
          XY, n = order, min_iterations = 1 + q, max_iterations = 5 + q,
          beta = .1, phi_F_XY = .5, q_F_XY = q, tails_count_threshold = threshold)
        fit <- ddfs(XY, n = order, min.intknots = 1, max.intknots = 5,
                    beta = .1, phi_F = .5, q_F = q, tails_count_threshold = threshold)
        expect_s3_class(fit, "ddfs")
        expect_identical(names(fit), c(names(legacy), "extcall"))
        expect_equal(fit[names(legacy)], legacy, tolerance = 1e-9)
        class(legacy) <- "ddfs"
        query <- XY[1:10, , drop = FALSE]
        expect_equal(d.ddfs(fit, query), d.ddfs(legacy, query), tolerance = 1e-9)
        expect_equal(p.ddfs(fit, query), p.ddfs(legacy, query), tolerance = 1e-9)
        expect_equal(coef(fit), coef(legacy), tolerance = 1e-9)
        expect_equal(knots(fit, options = "internal"),
                     knots(legacy, options = "internal"), tolerance = 1e-10)
      }
    }
  }
})

test_that("bivariate dispatch uses the new engine and preserves automatic tuning", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  set.seed(932)
  XY <- data.frame(first = runif(100), second = runif(100))
  testthat::local_mocked_bindings(
    BivariateDensityFitter = function(...) stop("Legacy fitter must not be called"),
    choose_params = function(...) list(min.intknots = 1, beta = .1, phi_F = .95, q_F = 1),
    .package = "ddfs")
  expect_message(fit <- ddfs(XY, n = 3, max.intknots = 2, tails_count_threshold = 0),
                 "Automatic parameter selection applied")
  expect_identical(fit$type, "Biv - DDFS")
  expect_identical(fit$args$XY, XY)
  expect_output(print(fit), "Biv - DDFS")
  expect_output(summary(fit), "Biv - DDFS")
  expect_equal(predict(fit), as.numeric(fit$f_XY_hat$pred))
  expect_error(ddfs(XY, n = 3, min.intknots = 1, beta = .1, phi_F = .95, q_F = 1,
                    max.coef = 15), "Tensor basis requires")
})

test_that("a stopping threshold at the iteration cap retains the final fit", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  set.seed(933)
  XY <- cbind(runif(100), runif(100))
  legacy <- ddfs:::BivariateDensityFitter(
    XY, n = 3, min_iterations = 3, max_iterations = 4,
    phi_F_XY = 1e-10, q_F_XY = 2, tails_count_threshold = 0)
  fit <- ddfs(XY, n = 3, min.intknots = 1, max.intknots = 2,
              beta = .1, phi_F = 1e-10, q_F = 2, tails_count_threshold = 0)
  expect_equal(fit$model, 4)
  expect_equal(fit[names(legacy)], legacy, tolerance = 1e-9)
})

test_that("bivariate plots and optional residual displays remain available", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  set.seed(934)
  XY <- cbind(runif(100), runif(100))
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  on.exit({ grDevices::dev.off(); unlink(output) }, add = TRUE)
  expect_no_error(fit <- ddfs(
    XY, n = 3, min.intknots = 1, max.intknots = 1,
    beta = .1, phi_F = .95, q_F = 1, tails_count_threshold = 0,
    plot = TRUE, resids_plot = TRUE, stop_plotting = 2))
  expect_no_error(plot(fit, fit = "pdf"))
  expect_no_error(plot(fit, fit = "cdf"))
})
