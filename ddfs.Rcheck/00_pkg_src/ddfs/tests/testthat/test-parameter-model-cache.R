test_that("selection models are reused and changed files are reloaded", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  testthat::local_mocked_bindings(
    .ddfs_parameter_models = new.env(parent = emptyenv()), .package = "ddfs")
  # Environment identity distinguishes reuse from a second deserialization.
  model <- new.env(parent = emptyenv())
  model$value <- 1
  saveRDS(model, path)
  first <- ddfs:::read_parameter_model(path)
  expect_identical(ddfs:::read_parameter_model(path), first)
  # Changing serialized size invalidates the cache even on coarse file clocks.
  saveRDS(list(value = seq_len(1000)), path, compress = FALSE)
  expect_identical(ddfs:::read_parameter_model(path), readRDS(path))
  unlink(path)
  expect_error(suppressWarnings(ddfs:::read_parameter_model(path)), "cannot open")
  saveRDS(list(value = 3), path)
  expect_identical(ddfs:::read_parameter_model(path), list(value = 3))
})

test_that("different model files and returned lists remain independent", {
  paths <- c(tempfile(fileext = ".rds"), tempfile(fileext = ".rds"))
  on.exit(unlink(paths), add = TRUE)
  testthat::local_mocked_bindings(
    .ddfs_parameter_models = new.env(parent = emptyenv()), .package = "ddfs")
  saveRDS(list(value = 1), paths[1])
  saveRDS(list(value = 2), paths[2])
  first <- ddfs:::read_parameter_model(paths[1])
  first$value <- 99
  expect_identical(ddfs:::read_parameter_model(paths[1]), list(value = 1))
  expect_identical(ddfs:::read_parameter_model(paths[2]), list(value = 2))
})

test_that("cold and warm caches preserve automatic selection exactly", {
  testthat::local_mocked_bindings(
    .ddfs_parameter_models = new.env(parent = emptyenv()), .package = "ddfs")
  set.seed(960)
  for (size in c(100, 600)) {
    for (dimension in 1:2) {
      X <- if (dimension == 1) rnorm(size) else matrix(rnorm(2 * size), ncol = 2)
      for (type in (if (dimension == 1) c("bw_top_mean_Q3", "default_top_mean_Q3") else "bw_top_mean_Q3")) {
        reference <- local({
          testthat::local_mocked_bindings(read_parameter_model = base::readRDS, .package = "ddfs")
          suppressWarnings(ddfs:::choose_params(X, type))
        })
        expect_identical(suppressWarnings(ddfs:::choose_params(X, type)), reference)
        expect_identical(suppressWarnings(ddfs:::choose_params(X, type)), reference)
      }
    }
  }
})

test_that("cached automatic selection preserves complete fitted objects", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  testthat::local_mocked_bindings(
    .ddfs_parameter_models = new.env(parent = emptyenv()), .package = "ddfs")
  set.seed(961)
  for (X in list(rnorm(100), rexp(100))) {
    fit_once <- function() suppressMessages(ddfs(X, n = 3))
    reference <- local({
      testthat::local_mocked_bindings(read_parameter_model = base::readRDS, .package = "ddfs")
      fit_once()
    })
    expect_identical(fit_once(), reference)
    expect_identical(fit_once(), reference)
  }
})

test_that("cached bivariate automatic fits agree exactly", {
  skip_if_not(exists("stageALoopND", asNamespace("GeDS"), inherits = FALSE))
  for (size in c(100, 1000)) {
    for (parameter in c("bw", "q_F", "phi_F_1", "phi_F_2", "beta")) {
      path <- paste0("choose_", parameter, "_N=", size, "_bw_top_mean_Q3_biv.rds")
      expect_true(nzchar(system.file(path, package = "ddfs")), info = path)
    }
  }
  testthat::local_mocked_bindings(
    .ddfs_parameter_models = new.env(parent = emptyenv()), .package = "ddfs")
  set.seed(962)
  X <- cbind(rnorm(100), rnorm(100))
  fit_once <- function() suppressMessages(ddfs(X, n = 3))
  reference <- local({
    testthat::local_mocked_bindings(read_parameter_model = base::readRDS, .package = "ddfs")
    fit_once()
  })
  expect_identical(fit_once(), reference)
  expect_identical(fit_once(), reference)
})
