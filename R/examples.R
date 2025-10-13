#' @title Simulate Data from Benchmark Distributions
#'
#' @description
#' Simulates univariate data from a variety of standard and benchmark distributions
#' used in nonparametric density estimation, such as normal, Student-\emph{t},
#' mixtures of Gaussians, heavy-tailed, skewed, or multimodal distributions.
#'
#' Many of the examples were extracted from in Cui et al. (2020) and
#' other literature, and is primarily intended for simulation studies and
#' performance comparisons in density estimation.
#'
#' @param N Integer. Number of observations to simulate.
#' @param ex Character string. Specifies the example/distribution to simulate from.
#' Supported options include:
#' \itemize{
#'   \item \code{"Gaussian"}: Standard normal distribution \eqn{N(0,1)}
#'   \item \code{"Student-t"}: Student-\eqn{t} distribution with 6 degrees of freedom
#'   \item \code{"Exponential"}: \eqn{\text{Exp}(\lambda=1)}, support \eqn{x \geq 0}
#'   \item \code{"Chi-square"}: \eqn{\chi^2_4}, support \eqn{x > 0}
#'   \item \code{"Gamma"}: \eqn{\text{Gamma}(k=9, \theta=0.5)}, support \eqn{x \geq 0}
#'   \item \code{"Weibull"}: \eqn{\text{Weibull}(\lambda=1, k=5)}, support \eqn{x \geq 0}
#'   \item \code{"Log-normal"}: \eqn{\text{Lognormal}(\mu=0, \sigma=1)}, support \eqn{x > 0}
#'   \item \code{"Nakagami"}: Nakagami distribution with \eqn{\mu = \omega = 2}
#'
#'   \item \code{"Merton's jump diffusion"}: Merton model with \eqn{\sigma = 0.08}, \eqn{\lambda = 3}, \eqn{\mu_J = -0.01}, \eqn{\sigma_J = 0.4}, \eqn{\Delta t = 1/4}
#'   \item \code{"Kou's double exponential"}: Kou model with \eqn{\sigma = 0.04}, \eqn{\lambda = 2}, \eqn{p_{\text{up}} = 0.4}, \eqn{\eta_1 = 3}, \eqn{\eta_2 = 5}, \eqn{\Delta t = 1/4}
#'
#'   \item \code{"GEV Type I"}: Generalized Extreme Value (Gumbel) with shape \eqn{\xi = 0}
#'   \item \code{"GEV Type II"}: Generalized Extreme Value (Fréchet) with \eqn{\xi = 0.5}
#'   \item \code{"GEV Type III"}: Generalized Extreme Value (Weibull) with \eqn{\xi = -0.5}
#'
#'   \item \code{"Kurtotic unimodal"}: \eqn{\frac{2}{3}N(0,1) + \frac{1}{3}N(0, (1/10)^2)}
#'   \item \code{"Outlier"}: \eqn{\frac{1}{10}N(0,1) + \frac{9}{10}N(0, (1/10)^2)}
#'   \item \code{"Skewed unimodal"}: \eqn{\frac{1}{5}N(0,1) + \frac{1}{5}N(0.5, (2/3)^2) + \frac{3}{5}N(13/12, (5/9)^2)}
#'   \item \code{"Strongly skewed"}: \eqn{\sum_{k=0}^7 \frac{1}{8}N(3((2/3)^k - 1), (2/3)^{2k})}
#'   \item \code{"Bimodal"}: \eqn{\frac{1}{2}N(0, (1/10)^2) + \frac{1}{2}N(5, 1)}
#'   \item \code{"Separated bimodal"}: \eqn{\frac{1}{2}N(-2, (1/2)^2) + \frac{1}{2}N(2, (1/2)^2)}
#'   \item \code{"Skewed bimodal"}: \eqn{\frac{3}{4}N(0,1) + \frac{1}{4}N(3/2, (1/3)^2)}
#'   \item \code{"MixGauss"}: \eqn{0.15N(-0.25,1/3) + 0.85N(3.25,1)}
#'   \item \code{"MixGauss2"}: \eqn{\frac{5}{6}N(3,1) + \frac{5}{36}N(8, (1/3)^2) + \frac{1}{36}N(10, (1/9)^2)}
#'   \item \code{"MixGauss3"}: \eqn{0.3N(1, 0.5^2) + 0.2N(4.5, 1.2^2) + 0.5N(8, 0.8^2)}
#'   \item \code{"Mix1d"}: \eqn{0.8\chi^2(3) + 0.2N(7,1)}
#'   \item \code{"Trimodal"}: \eqn{\frac{1}{3}\sum_{k=0}^2 N(80k, (k+1)^4)}
#'   \item \code{"Smooth comb"}: \eqn{\sum_{k=0}^5 \frac{2^{5-k}}{63}N\left(\frac{65 - 96/2^k}{21}, \left(\frac{32/63}{2^k}\right)^2\right)}
#'   \item \code{"Claw"}: \eqn{\frac{1}{2}N(0,1) + \sum_{k=0}^4 \frac{1}{10}N(k/2 - 1, (1/10)^2)}
#'
#'   \item \code{"BivGauss_0"}: Bivariate standard normal distribution with correlation \eqn{\rho = 0}
#'   \item \code{"BivGauss_0.2"}: Bivariate standard normal distribution with correlation \eqn{\rho = 0.2}
#'   \item \code{"BivGauss_0.5"}: Bivariate standard normal distribution with correlation \eqn{\rho = 0.5}
#'   \item \code{"BivGauss_0.8"}: Bivariate standard normal distribution with correlation \eqn{\rho = 0.8}
#'
#'   \item \code{"BivBiModGauss"}: Bivariate Gaussian mixture: \eqn{0.8 \cdot N(\mu_1, \Sigma_1) + 0.2 \cdot N(\mu_2, \Sigma_2)}, with
#'   \eqn{\mu_1 = (5.5, 5.5)}, \eqn{\Sigma_1 = \begin{pmatrix} 0.36 & 0.108 \\ 0.108 & 0.36 \end{pmatrix}},
#'   \eqn{\mu_2 = (7, 7)}, and \eqn{\Sigma_2 = \begin{pmatrix} 0.16 & 0.048 \\ 0.048 & 0.16 \end{pmatrix}}
#'
#'   \item \code{"BivGaussSkewed"}: \eqn{\frac{1}{5} N(0,0,1,1,0) + \frac{1}{5} N(\frac{1}{2},\frac{1}{2},(\frac{2}{3})^2,(\frac{2}{3})^2,0) + \frac{3}{5} N(\frac{13}{12},\frac{13}{12},(\frac{5}{9})^2,(\frac{5}{9})^2,0)}

#' }
#'
#' @return A named \code{list} with the following elements:
#' \describe{
#'   \item{\code{X}}{The sorted vector of simulated observations.}
#'   \item{\code{f_X}}{The true density evaluated at the simulated \code{X}.}
#'   \item{\code{F_X}}{The true cumulative distribution function evaluated at \code{X}.}
#'   \item{\code{f_X_func}}{A function \code{f(x)} to evaluate the true pdf at arbitrary points.}
#'   \item{\code{F_X_func}}{A function \code{F(x)} to evaluate the true cdf at arbitrary points.}
#' }
#'
#' @details
#' All generated distributions are continuous and generated using inverse
#' transform sampling is used over the appropriate quantile range.
#'
#' @examples
#' N <- 10000
#' examples <- c("Gaussian", "Student-t", "Exponential", "Chi-square", "Gamma", "Weibull",
#' "Log-normal", "Nakagami", "Kurtotic unimodal", "Outlier", "Skewed unimodal",
#' "Strongly skewed", "Merton's jump diffusion", "Kou's double exponential",
#' "GEV Type I", "GEV Type II", "GEV Type III",
#' "MixGauss", "Mix1d", "MixGauss2",
#' "Bimodal", "Separated bimodal", "Skewed bimodal",
#' "Trimodal", "Smooth comb", "Claw")
#' par(mfrow = c(2,3))
#' for (ex in examples) {
#'
#'   set.seed(123)
#'   sim <- sim.dist(N, ex = ex)
#'
#'   # pdf
#'   f_X <- sim$f_X
#'   hist(sim$X, breaks = 100, col = rgb(0.2, 0.4, 0.6, 0.6), border = "black",
#'        main = ex,
#'        xlab = "X values", ylab = "Frequency",
#'        freq = FALSE)  # Set freq = FALSE for density scale
#'   lines(sim$X, f_X, col = "darkblue", lwd = 2)
#'
#'   # cdf
#'   F_X <- sim$F_X
#'   plot(sim$X, cumsum(rep(1/N, N)))
#'   lines(sim$X, F_X)
#' }
#'
#' @import stats
#' @importFrom extRemes devd pevd revd
#' @importFrom mvtnorm dmvnorm pmvnorm rmvnorm
#' @importFrom nakagami dnaka pnaka rnaka
#'
#' @references
#' Cui, Z., Kirkby, J. L., & Nguyen, D. (2020). Nonparametric density estimation by B-spline duality.
#' \emph{Econometric Theory}, 36(2), 250–291. \doi{10.1017/S0266466619000112}

#' @export

sim.dist <- function(N, ex) {

  if (!(ex %in% c("Gaussian", "Student-t", "Exponential", "Chi-square", "Gamma", "Weibull",
                  "Log-normal", "Nakagami", "Kurtotic unimodal", "Outlier", "Skewed unimodal",
                  "Strongly skewed", "Merton's jump diffusion", "Kou's double exponential",
                  "GEV Type I", "GEV Type II", "GEV Type III",
                  "MixGauss", "Mix1d", "MixGauss2", "MixGauss3",
                  "Bimodal", "Separated bimodal", "Skewed bimodal",
                  "Trimodal", "Smooth comb", "Claw",
                  "BivGauss_0", "BivGauss_0.2", "BivGauss_0.5", "BivGauss_0.8", "BivBiModGauss",
                  "BivGaussSkewed", "BivGaussKurtotic",
                  "BivBiModGauss2", "BivBiModGauss3", "BivBiModGauss4", "BivBiModGauss5",
                  "BivTriModGauss1", "BivTriModGauss2", "BivTriModGauss3",
                  "BivQuadriModGauss"))) {
    stop("Invalid distribution name.")
  }

  if (ex == "Gaussian") {
    # N(0,1)
    # Parameters
    mu = 0; sigma = 1; a = -Inf; b = +Inf
    # Generate random samples for each component
    p_lower <- pnorm(a, mean = mu, sd = sigma)
    p_upper <- pnorm(b, mean = mu, sd = sigma)
    U <- runif(N, min = p_lower, max = p_upper)
    X <- sort(qnorm(U, mean = mu, sd = sigma))
    # pdf and cdf
    f_X_func <- function(x) dnorm(x, mean = mu, sd = sigma)
    f_X <- dnorm(X, mean = mu, sd = sigma)
    F_X_func <- function(x) pnorm(x, mean = mu, sd = sigma)
    F_X <- pnorm(X, mean = mu, sd = sigma)

  } else if (ex == "Student-t") {
    # t_v, v=6
    # Parameters
    df = 6; a = -Inf; b = +Inf
    # Generate random samples for each component
    p_lower <- pt(a, df = df)
    p_upper <- pt(b, df = df)
    U <- runif(N, min = p_lower, max = p_upper)
    X <- sort(qt(U, df = df))
    # pdf and cdf
    f_X_func <- function(x) dt(x, df = df)
    f_X <- dt(X, df = df)
    F_X_func <- function(x) pt(x, df = df)
    F_X <- pt(X, df = df)

  } else if (ex == "Exponential") {
    # λ = 1, x ≥ 0
    # Parameters
    lambda = 1; a = 0; b = +Inf
    # Generate random samples
    p_lower <- pexp(a, rate = lambda)
    p_upper <- pexp(b, rate = lambda)
    U <- runif(N, min = p_lower, max = p_upper)
    X <- sort(qexp(U, rate = lambda))
    # pdf and cdf
    f_X_func <- function(x) dexp(x, rate = lambda)
    f_X <- dexp(X, rate = lambda)
    F_X_func <- function(x) dexp(x, rate = lambda)
    F_X <- pexp(X, rate = lambda)

  } else if (ex == "Chi-square") {
    # χ2(k), k = 4, x > 0
    # Parameters
    df = 4; a = 0; b = +Inf
    # Generate random samples for each component
    p_lower <- pchisq(a, df = df)
    p_upper <- pchisq(b, df = df)
    U <- runif(N, min = p_lower, max = p_upper)
    X <- sort(qchisq(U, df = df))
    # pdf and cdf
    f_X_func <- function(x) dchisq(x, df = df)
    f_X <- dchisq(X, df = df)
    F_X_func <- function(x) pchisq(x, df = df)
    F_X <- pchisq(X, df = df)

  } else if (ex == "Gamma") {
    # x ≥ 0, k = 9 (shape), θ = 0.5 (scale)
    # Parameters
    k = 9; theta = 0.5; a = 0; b = +Inf
    # Generate random samples for each component
    p_lower <- pgamma(a, shape = k, scale = theta)
    p_upper <- pgamma(b, shape = k, scale = theta)
    U <- runif(N, min = p_lower, max = p_upper)
    X <- sort(qgamma(U, shape = k, scale = theta))
    # pdf and cdf
    f_X_func <- function(x) dgamma(x, shape = k, scale = theta)
    f_X <- dgamma(X, shape = k, scale = theta)
    F_X_func <- function(x) pgamma(x, shape = k, scale = theta)
    F_X <- pgamma(X, shape = k, scale = theta)

  } else if (ex == "Weibull") {
    # Parameters
    lambda = 1; k = 5; a = 0; b = +Inf
    # Generate random samples for each component
    p_lower <- pweibull(a, shape = k, scale = lambda)
    p_upper <- pweibull(b, shape = k, scale = lambda)
    U <- runif(N, min = p_lower, max = p_upper)
    X <- sort(qweibull(U, shape = k, scale = lambda))
    # pdf and cdf
    f_X_func <- function(x) dweibull(x, shape = k, scale = lambda)
    f_X <- dweibull(X, shape = k, scale = lambda)
    F_X_func <- function(x) pweibull(x, shape = k, scale = lambda)
    F_X <- pweibull(X, shape = k, scale = lambda)

  } else if (ex == "Log-normal") {
    # Parameters
    mu = 0; sigma = 1; a = 0; b = +Inf
    # Generate random samples for each component
    p_lower <- plnorm(a, meanlog = mu, sdlog = sigma)
    p_upper <- plnorm(b, meanlog = mu, sdlog = sigma)
    U <- runif(N, min = p_lower, max = p_upper)
    X <- sort(qlnorm(U, meanlog = mu, sdlog = sigma))
    # pdf and cdf
    f_X_func <- function(x) dlnorm(x, meanlog = mu, sdlog = sigma)
    f_X <- dlnorm(X, meanlog = mu, sdlog = sigma)
    F_X_func <- function(x) plnorm(x, meanlog = mu, sdlog = sigma)
    F_X <- plnorm(X)

  } else if (ex == "Nakagami") {
    mu = omega = 2
    X <- nakagami::rnaka(N, shape = mu, scale = omega)
    X <- sort(X)

    # pdf and cdf
    f_X_func <- function(x) dnaka(x, shape = mu, scale = omega)
    f_X <- nakagami::dnaka(X, shape = mu, scale = omega)
    F_X_func <- function(x) pnaka(x, shape = mu, scale = omega)
    F_X <- nakagami::pnaka(X, shape = mu, scale = omega)

  } else if (ex == "Merton's jump diffusion") {
    merton <- merton_sim(N_sim = N, sigma = 0.08, lambda  = 3, mu_J = -0.01,
                         sigma_J = 0.4, Delta_t = 1/4)
    X <- sort(merton$X)
    # pdf and cdf
    f_X_func <- merton$f_X_func
    f_X <- f_X_func(X)
    F_X_func <- merton$F_X_func
    F_X <- F_X_func(X)

  } else if (ex == "Kou's double exponential") {
    kou <- kou_sim(N_sim = N, Delta_t = 1/4,
                   sigma = 0.04, lambda  = 2, p_up = 0.4,
                   eta_1 = 3, eta_2 = 5)

    X <- sort(kou$X)
    # pdf and cdf
    f_X_func <- kou$f_X_func
    f_X <- f_X_func(X)
    F_X_func <- kou$F_X_func
    F_X <- F_X_func(X)

  } else if (ex == "GEV Type I" || ex == "GEV Type II" || ex == "GEV Type III") {

    mu = 0; sigma = 1
    if (ex == "GEV Type I") {
      ## Type I (Gumbel) (ξ = 0)
      xi = 0
    } else if (ex == "GEV Type II") {
      ## Type II (Fréchet) (ξ > 0)
      xi = 0.5
    } else if (ex == "GEV Type III") {
      ## Type III (Weibull) (ξ < 0)
      xi = -0.5
    }
    X <- sort(extRemes::revd(N, loc = mu , scale = sigma, shape = xi))

    # pdf and cdf
    f_X_func <- function(x) devd(x, mu, sigma, xi)
    f_X <- extRemes::devd(X, mu, sigma, xi)
    F_X_func <- function(x) pevd(x, mu, sigma, xi)
    F_X <- extRemes::pevd(X, loc = mu, scale = sigma, shape = xi)

  } else if (ex == "Kurtotic unimodal" || ex == "Skewed unimodal" ||
             ex == "Strongly skewed" || ex == "Outlier" || ex == "Smooth comb" ||
             ex == "Bimodal" || ex == "Separated bimodal" || ex == "Skewed bimodal" ||
             ex == "Trimodal" || ex == "Claw" || ex == "MixGauss" || ex == "MixGauss2" ||
             ex == "MixGauss3") {

    # Parameters
    a <- -Inf; b <- +Inf

    if (ex == "Kurtotic unimodal") {
      # (2/3)N(0,1) + (1/3)N(0,(1/10)^2)
      mu <- c(0, 0)
      sigma <- c(1, 1/10)
      p <- c(2/3, 1/3)

    } else if (ex == "Skewed unimodal") {
      # (1/5)N(0,1) + (1/5)N(1/2,(2/3)^2) + (3/5)N(13/12,(5/9)^2)
      mu <- c(1, 1/2, 13/12)
      sigma <- c(1, 2/3, 5/9)
      p <- c(1/5, 1/5, 3/5)

    } else if (ex == "Strongly skewed") {
      # \sum_{k=0}^7(1/8)N (3((2/3)^k − 1),(2/3)^2k)
      K <- 8
      k_vals <- 0:(K - 1)
      mu <- 3 * ((2/3)^k_vals - 1)
      sigma <- (2/3)^k_vals
      p <- rep(1 / K, K)

    } else if (ex == "Outlier") {
      # (1/10)N(0,1) + (9/10)N(0,(1/10)^2)
      mu <- c(0, 0)
      sigma <- c(1, 1/10)
      p <- c(1/10, 9/10)

    } else if (ex == "Smooth comb") {
      # \sum_{k=0}^5 2^{5-k}/63N((65-96/2^k)/21,((32/63)/2^k)^2)
      K <- 6
      k_vals <- 0:(K - 1)
      mu <- (65-96/2^k_vals)/21
      sigma <- (32/63)/2^k_vals
      p <- 2^(5-k_vals)/63

    } else if (ex == "Bimodal") {
      # (1/2)N(0,(1/10)^2) + (1/2)N(5,1)
      mu <- c(0, 5)
      sigma <- c(1/10, 1)
      p <- c(1/2, 1/2)

    } else if (ex == "Separated bimodal") {
      # (1/2)N(-2,(1/2)^2) + (1/2)N(2,(1/2)^2)
      mu <- c(-2, 2)
      sigma <- c(1/2, 1/2)
      p <- c(1/2, 1/2)

    } else if (ex == "Skewed bimodal") {
      # (3/4)N(0,1) + (1/4)N(3/2,(1/3)^2)
      mu <- c(0, 3/2)
      sigma <- c(1, 1/3)
      p <- c(3/4, 1/4)

    } else if (ex == "Trimodal") {
      # (1/3)\sum_{k=0}^2N(80k, (k+1)^4)
      K <- 3
      k_vals <- 0:(K - 1)
      mu <- 80*k_vals
      sigma <- (k_vals + 1)^2
      p <- rep(1/3, K)

    } else if (ex == "Claw") {
      # (1/2)N(0,1)+\sum_{k=0}^4(1/10)N(\frac{k}{2}-1, (\frac{1}{10})^2)
      K <- 5
      k_vals <- 0:(K - 1)
      mu <- c(0, k_vals/2-1)
      sigma <- c(1, rep(1/10, K) )
      p <- c(1/2, rep(1/10, K))

    } else if (ex == "MixGauss") {
      # 0.15N(−0.25, 1/3) + 0.85N(3.25, 1)
      mu <- c(-0.25, 3.25)
      sigma <- c(sqrt(1/3), sqrt(1))
      p <- c(0.15, 0.85)

    } else if (ex == "MixGauss2") {

      # (5/6)N(3, 1) + (5/36)N(8, (1/3)^2) + (1/36)N(10, (1/9)^2)
      mu <- c(3, 8, 10)
      sigma <- c(1, 1/3, 1/9)
      p <- c(5/6, 5/36, 1/36)

    } else if (ex == "MixGauss3") {
      # 0.3N(1, 0.5^2) + 0.2N(4.5, 1.2^2) + 0.5N(8, 0.8^2)
      mu <- c(1, 4.5, 8)
      sigma <- c(sqrt(0.5^2), sqrt(1.2^2), sqrt(0.8^2))
      p <- c(0.3, 0.2, 0.5)
    }


    # Generate random samples for each component
    X <- vector("list", length(p))
    for (i in 1:length(p)) {
      p_lower <- pnorm(a, mean = mu[i], sd = sigma[i])
      p_upper <- pnorm(b, mean = mu[i], sd = sigma[i])
      U <- runif(N, min = p_lower, max = p_upper)
      X[[i]] <- qnorm(U, mean = mu[i], sd = sigma[i])
    }
    names(X) <- paste0("X", 1:length(p))
    # Generate the mixture components
    z <- t(rmultinom(N, 1, p))
    # Combine the components
    X <- rowSums(z * do.call(cbind, X))
    X <- sort(X)

    # pdf
    f_X_func <- (function(p, mu, sigma) {
      function(x) rowSums(
        vapply(1:length(p), function(k) {
          p[k] * dnorm(x, mean = mu[k], sd = sigma[k])
        }, numeric(length(x)))
      )
    })(p, mu, sigma)
    f_X <- rowSums(
      vapply(1:length(p), function(k) {
        p[k] * dnorm(X, mean = mu[k], sd = sigma[k])
      }, numeric(length(X)))
    )
    # cdf
    F_X_func <- (function(p, mu, sigma) {
      function(x) rowSums(
        vapply(1:length(p), function(k) {
          p[k] * pnorm(x, mean = mu[k], sd = sigma[k])
        }, numeric(length(x)))
      )
    })(p, mu, sigma)
    F_X <- rowSums(
      vapply(1:length(p), function(k) {
        p[k] * pnorm(X, mean = mu[k], sd = sigma[k])
      }, numeric(length(X)))
    )


  } else if (ex == "Mix1d") {

    # 0.8χ2(3) + 0.2N(7, 1) in [0,10]
    df = 3; p = 0.8
    mu2 = 7; sigma2 = 1
    a1 = 0; b1 = +Inf
    a2 = -Inf; b2 = +Inf
    p_lower1 <- pchisq(a1, df = df)
    p_upper1 <- pchisq(b1, df = df)
    U1 <- runif(N, min = p_lower1, max = p_upper1)
    p_lower2 <- pnorm(a2, mean = mu2, sd = sigma2)
    p_upper2 <- pnorm(b2, mean = mu2, sd = sigma2)
    U2 <- runif(N, min = p_lower2, max = p_upper2)

    X1 <- sort(qchisq(U1, df = df))
    X2 <- sort(qnorm(U2, mean = mu2, sd = sigma2))

    z <- rbinom(N, size = 1, prob = p)
    X <- z * X1 + (1 - z) * X2
    X <- sort(X)

    # pdf and cdf
    f_X_func <- function(x) p*dchisq(x, df = df) + (1-p)*dnorm(x, mu2, sigma2)
    f_X <- p*dchisq(X, df = df) + (1-p)*dnorm(X, mu2, sigma2)
    F_X_func <- function(x) p*pchisq(x, df = df) + (1-p)*pnorm(x, mu2, sigma2)
    F_X <- p*pchisq(X, df = df) + (1-p)*pnorm(X, mu2, sigma2)

  } else if (ex == "BivGauss_0" || ex == "BivGauss_0.2" || ex == "BivGauss_0.5" || ex == "BivGauss_0.8") {

    # Bivariate Gaussian with correlation rho
    if (ex == "BivGauss_0") {
      rho <- 0
    } else if (ex == "BivGauss_0.2") {
      rho <- 0.2
    } else if (ex == "BivGauss_0.5") {
      rho <- 0.5
    } else if (ex == "BivGauss_0.8") {
      rho <- 0.8
    }

    Sigma <- matrix(c(1, rho, rho, 1), 2)
    X <- rmvnorm(N, mean = c(0,0), sigma = Sigma)

    f_X_func <- function(x) dmvnorm(x, mean = c(0,0), sigma = Sigma)
    f_X <- dmvnorm(X, mean = c(0,0), sigma = Sigma)

    F_X_func <- function(x) apply(x, 1, function(row)
      pmvnorm(upper = row, mean = c(0,0), sigma = Sigma)[1]
    )
    F_X <- apply(X, 1, function(row)
      pmvnorm(upper = row, mean = c(0,0), sigma = Sigma)[1]
    )

  } else if (ex == "BivBiModGauss") {
    # Parameters
    a1 <- a2 <- 0; b1 <- b2 <- 9
    p <- 0.8
    # Define parameters for the first component
    mu1 <- c(5.5, 5.5)
    sigma1_x <- 0.6; sigma1_y <- 0.6; rho1 <- 0.3
    sigma1 <- matrix(c(sigma1_x^2, rho1*sigma1_x*sigma1_y, rho1*sigma1_y*sigma1_x, sigma1_y^2), ncol=2)
    # Define parameters for the second component
    mu2 <- c(7, 7)
    sigma2_x <- 0.4; sigma2_y <- 0.4; rho2 <- 0.3
    sigma2 <- matrix(c(sigma2_x^2, rho2*sigma2_x*sigma2_y, rho2*sigma2_y*sigma2_x, sigma2_y^2), ncol=2)

    # Compute K (the normalization constant)
    K1 <- pmvnorm(lower = c(a1, a2), upper = c(b1, b2),
                  mean = mu1,
                  sigma = sigma1)[1]
    K2 <- pmvnorm(lower = c(a1, a2), upper = c(b1, b2),
                  mean = mu2,
                  sigma = sigma2)[1]
    K <- p * K1 + (1-p)* K2

    # Simulate data
    XY1 <- rmvnorm(n=N, mean=mu1, sigma=sigma1)
    XY2 <- rmvnorm(n=N, mean=mu2, sigma=sigma2)
    z <- rbinom(N, size = 1, prob = p)
    X <- z * XY1 + (1 - z) * XY2

    f_X_func <- function(x) {
      g_XY <- p * mapply(function(xi, yi) dmvnorm(c(xi, yi), mean = mu1, sigma = sigma1), x[,1], x[,2]) +
        (1-p) * mapply(function(xi, yi) dmvnorm(c(xi, yi), mean = mu2, sigma = sigma2), x[,1], x[,2])
      return(g_XY/K)
    }
    f_X <- f_X_func(X)

    F_X_func <- function(x) {
      G_XY <- p * mapply(function(xi, yi)
        pmvnorm(upper = c(xi, yi), mean = mu1, sigma = sigma1),
        x[,1], x[,2]) +
        (1 - p) * mapply(function(xi, yi)
          pmvnorm(upper = c(xi, yi), mean = mu2, sigma = sigma2),
          x[,1], x[,2])
      return(G_XY/K)
    }

    F_X <- F_X_func(X)

  } else if (ex %in% c("BivGaussSkewed", "BivGaussKurtotic",
                       "BivBiModGauss2", "BivBiModGauss3", "BivBiModGauss4", "BivBiModGauss5",
                       "BivTriModGauss1", "BivTriModGauss2", "BivTriModGauss3",
                       "BivQuadriModGauss") ){

    # Parameters
    if (ex == "BivGaussSkewed") {
      mu    <- c(0, 0, 1/2, 1/2, 13/12, 13/12)  # flattened means
      sigma <- c(1, 1, 2/3, 2/3, 5/9, 5/9)      # flattened sds
      rho   <- c(0, 0, 0)                       # correlations
      p     <- c(1/5, 1/5, 3/5)                 # mixture weights
    } else if (ex == "BivGaussKurtotic") {
      mu    <- c(0, 0, 0, 0)                    # flattened means
      sigma <- c(1, 2, 2/3, 1/3)                # flattened sds
      rho   <- c(1/2, -1/2)                     # correlations
      p     <- c(2/3, 1/3)                      # mixture weights

    } else if (ex == "BivBiModGauss2") {
      mu    <- c(-1, 0, 1, 0)                  # flattened means
      sigma <- c(2/3, 2/3, 2/3, 2/3)           # flattened sds
      rho   <- c(0, 0)                         # correlations
      p     <- c(1/2, 1/2)                     # mixture weights

    } else if (ex == "BivBiModGauss3") {
      mu    <- c(-3/2, 0, 3/2, 0)             # flattened means
      sigma <- c(1/4, 1, 1/4, 1)              # flattened sds
      rho   <- c(0, 0)                        # correlations
      p     <- c(1/2, 1/2)                    # mixture weights

    } else if (ex == "BivBiModGauss4") {
      mu    <- c(-1, 1, 1, -1)  # flattened means
      sigma <- c(2/3, 2/3, 2/3, 2/3)          # flattened sds
      rho   <- c(3/5, 3/5)                    # correlations
      p     <- c(1/2, 1/2)                    # mixture weights

    } else if (ex == "BivBiModGauss5") {
      mu    <- c(1, -1, -1, 1)  # flattened means
      sigma <- c(2/3, 2/3, 2/3, 2/3)          # flattened sds
      rho   <- c(7/10, 0)                     # correlations
      p     <- c(1/2, 1/2)                    # mixture weights

    } else if (ex == "BivTriModGauss1") {
      mu    <- c(-6/5, 6/5, 6/5, -6/5, 0, 0)      # flattened means
      sigma <- c(3/5, 3/5, 3/5, 3/5, 1/4, 1/4)    # flattened sds
      rho   <- c(3/10, -3/5, 1/5)                 # correlations
      p     <- c(9/20, 9/20, 1/10)                # mixture weights

    } else if (ex == "BivTriModGauss2") {
      mu    <- c(-6/5, 0, 6/5, 0, 0, 0)            # flattened means
      sigma <- c(3/5, 3/5, 3/5, 3/5, 3/5, 3/5)     # flattened sds
      rho   <- c(7/10, 7/10, -7/10)                # correlations
      p     <- c(1/3, 1/3, 1/3)                    # mixture weights

    } else if (ex == "BivTriModGauss3") {
      mu    <- c(-1, 0, 1, (2/3)*sqrt(3), 1, -(2/3)*sqrt(3))  # flattened means
      sigma <- c(3/5, 7/10, 3/5, 7/10, 3/5, 7/10)             # flattened sds
      rho   <- c(3/5, 0, 0)                                   # correlations
      p     <- c(3/7 , 3/7, 1/7)                              # mixture weights

    } else if (ex == "BivQuadriModGauss") {
      mu    <- c(-1, 1, -1, -1, 1, -1, 1, 1)                  # flattened means
      sigma <- c(2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3, 2/3)      # flattened sds
      rho   <- c(2/5, 3/5, -7/10, -1/2)                       # correlations
      p     <- c(1/8, 3/8, 1/8, 3/8)                          # mixture weights

    }

    # Reshape into list of components
    K <- length(p)
    mu_list <- split(matrix(mu, ncol = 2, byrow = TRUE), 1:K)
    sigma_list <- split(matrix(sigma, ncol = 2, byrow = TRUE), 1:K)

    Sigma_list <- lapply(1:K, function(k) {
      sd1 <- sigma_list[[k]][1]
      sd2 <- sigma_list[[k]][2]
      matrix(c(sd1^2, rho[k]*sd1*sd2,
               rho[k]*sd1*sd2, sd2^2), 2, 2)
    })

    ## --- Sampling from mixture ---------------------------------------
    rmvnorm_mixture <- function(N) {
      comp <- sample(1:K, size = N, replace = TRUE, prob = p)
      X <- matrix(NA_real_, nrow = N, ncol = 2)
      for (k in 1:K) {
        idx <- which(comp == k)
        if (length(idx) > 0) {
          X[idx, ] <- rmvnorm(length(idx), mean = mu_list[[k]], sigma = Sigma_list[[k]])
        }
      }
      X
    }

    ## --- Density of mixture ------------------------------------------
    f_X_func <- function(x) {
      if (is.vector(x)) x <- matrix(x, ncol = 2, byrow = TRUE)
      rowSums(sapply(1:K, function(k)
        p[k] * dmvnorm(x, mean = mu_list[[k]], sigma = Sigma_list[[k]])))
    }

    ## --- CDF of mixture ----------------------------------------------
    F_X_func <- function(x) {
      if (is.vector(x)) x <- matrix(x, ncol = 2, byrow = TRUE)
      apply(x, 1, function(row) {
        sum(sapply(1:K, function(k)
          p[k] * as.numeric(pmvnorm(upper = row,
                                    mean = mu_list[[k]],
                                    sigma = Sigma_list[[k]]))
        ))
      })
    }

    X <- rmvnorm_mixture(N)
    f_X <- f_X_func(X)
    F_X <- F_X_func(X)

  }


  return(list(X = X, f_X = f_X, F_X = F_X,
              f_X_func = f_X_func, F_X_func = F_X_func))

}

################################################################################
#################### Merton/Kou Jump Diffusion Process #########################
################################################################################

###############
## 1. Merton ##
###############
# X(Δt) = ln(S(t + Δt)/S(t)) = θΔt +σ W(Δt)+\sum^{N(Δt)}_{k=1} Jk; Jk ~ N(μ_J, σ_J)
# https://quant-next.com/the-merton-jump-diffusion-model/
merton_sim <- function(N_sim = 100, Delta_t = 1/4, type = "Merton",
                       sigma, lambda,
                       mu_J, sigma_J,
                       p_up, eta_1, eta_2) {

  # Check core parameters
  if (missing(sigma) || !is.numeric(sigma) || length(sigma) != 1) stop("sigma is missing, not numeric, or not scalar")
  if (missing(lambda) || !is.numeric(lambda) || length(lambda) != 1) stop("lambda is missing, not numeric, or not scalar")

  # Check type-specific parameters
  if (type == "Merton") {
    if (missing(mu_J) || !is.numeric(mu_J) || length(mu_J) != 1) stop("mu_J is missing, not numeric, or not scalar")
    if (missing(sigma_J) || !is.numeric(sigma_J) || length(sigma_J) != 1) stop("sigma_J is missing, not numeric, or not scalar")
  } else if (type == "Kou") {
    if (missing(p_up) || !is.numeric(p_up) || length(p_up) != 1) stop("p_up is missing, not numeric, or not scalar")
    if (missing(eta_1) || !is.numeric(eta_1) || length(eta_1) != 1) stop("eta_1 is missing, not numeric, or not scalar")
    if (missing(eta_2) || !is.numeric(eta_2) || length(eta_2) != 1) stop("eta_2 is missing, not numeric, or not scalar")
  } else {
    stop("type must be either 'Merton' or 'Kou'")
  }

  # f(x) = \sum_{n=0}^∞ P(N(Δt) = n) * f_N(x); N(Δt) ~ Poisson(λ*Δt)
  n_max <- qpois(0.99, lambda = lambda*Delta_t) # generate enough summands
  n <- 0:n_max
  p_n <- dpois(n, lambda = lambda*Delta_t)

  # Base drift parameter for the diffusion (often mu is the risk-free rate or total drift)
  mu <- 0

  # θ is chosen so that the process has the desired drift (this version ensures the martingale property)
  # k = E(Y_i-1)=e^{μ_J+0.5*σ_J^2}-1
  k <- exp(mu_J + 0.5 * sigma_J^2) - 1
  # θ = μ - 0.5*σ^2 - λk
  theta <- mu - 0.5 * sigma^2 - lambda * k
  # The (conditional) density of X[i] given the jump is Normal with
  # mean = θ*Δt + N(Δt)*μ_J and variance = σ^2*Δt + N(Δt)*σ_J^2.
  mu_cond <- theta * Delta_t + n * mu_J
  sigma_cond <- sqrt(sigma^2 * Delta_t + n * sigma_J^2)

  # Preallocate vectors for the simulated log-returns and their density values
  X <- numeric(N_sim)
  f_X <- numeric(N_sim)
  F_X <- numeric(N_sim)

  # 1. Simulate the diffusion (Brownian) component W(Δt) ~ √Δt·N(0,1) = √Δt·z
  z <- rnorm(N_sim)
  # 2. Simulate the number of jumps over the interval N(Δt)
  N_Delta_t <- rpois(N_sim, lambda * Delta_t)

  for (i in 1:N_sim) {

    # σ·W(Δt)
    diffusion <- sigma * sqrt(Delta_t) * z[i]

    # 3. Sum the jump sizes \sum_{k=1}^N(Δt) J_k
    sum_J_k <- sum(rnorm(N_Delta_t[i], mean = mu_J, sd = sigma_J))

    # 4. Combine the drift, diffusion, and jump components to get the log-return
    # θ·Δt + σ·W(Δt) + \sum_{k=1}^N(Δt) J_k
    X[i] <- theta * Delta_t + diffusion + sum_J_k

    # 5. Compute pdf and cdf of ln(Xt/X0)
    f_N <- dnorm(X[i], mean = mu_cond, sd = sigma_cond)
    f_X[i] <- sum(p_n*f_N)
    F_N <- pnorm(X[i], mean = mu_cond, sd = sigma_cond)
    F_X[i] <- sum(p_n*F_N)

  }

  # Save pdf and cdf functions
  f_X_func <- (function(p_n, mu_cond, sigma_cond) {
    function(x) {
      sapply(x, function(xi) {
        sum(p_n * dnorm(xi, mean = mu_cond, sd = sigma_cond))
      })
    }
  })(p_n, mu_cond, sigma_cond)

  F_X_func <- (function(p_n, mu_cond, sigma_cond) {
    function(x) {
      sapply(x, function(xi) {
        sum(p_n * pnorm(xi, mean = mu_cond, sd = sigma_cond))
      })
    }
  })(p_n, mu_cond, sigma_cond)

  return(list(X = X,
              f_X = f_X, f_X_func = f_X_func,
              F_X = F_X, F_X_func = F_X_func))
}

############
## 2. Kou ##
############
# X(Δt) = ln(S(t + Δt)/S(t)) = θΔt +σ W(Δt)+\sum^{N(Δt)}_{k=1} Jk;
# Jk ~ f_J(y) = p_{up}η1exp(−η1y)1[y ≥ 0] + (1 − p_{up})η2exp(η2y)1[y < 0]
# https://papers.ssrn.com/sol3/papers.cfm?abstract_id=242367
# https://www.sciencedirect.com/science/article/pii/S0898122108003477
kou_sim <- function(N_sim, Delta_t = 1/4,
                    sigma, lambda,
                    p_up, eta_1, eta_2) {

  # Check core parameters
  if (missing(sigma) || !is.numeric(sigma) || length(sigma) != 1) stop("sigma is missing, not numeric, or not scalar")
  if (missing(lambda) || !is.numeric(lambda) || length(lambda) != 1) stop("lambda is missing, not numeric, or not scalar")
  if (missing(p_up) || !is.numeric(p_up) || length(p_up) != 1) stop("p_up is missing, not numeric, or not scalar")
  if (missing(eta_1) || !is.numeric(eta_1) || length(eta_1) != 1) stop("eta_1 is missing, not numeric, or not scalar")
  if (missing(eta_2) || !is.numeric(eta_2) || length(eta_2) != 1) stop("eta_2 is missing, not numeric, or not scalar")

  # Base drift parameter for the diffusion (often mu is the risk-free rate or total drift)
  mu <- 0
  # θ = μ - 0.5*σ^2
  theta <- mu - 0.5 * sigma^2

  # Preallocate vectors for the simulated log-returns and their density values
  X <- numeric(N_sim)
  f_X <- numeric(N_sim)
  F_X <- numeric(N_sim)

  # 1. Simulate the diffusion (Brownian) component W(Δt) ~ √Δt·N(0,1) = √Δt·z
  z <- rnorm(N_sim)
  # 2. Simulate the number of jumps over the interval Delta_t
  N_Delta_t <- rpois(N_sim, lambda * Delta_t)

  for (i in 1:N_sim) {

    # σ·W(Δt)
    diffusion <- sigma * sqrt(Delta_t) * z[i]

    # 3. Sum the jump sizes \sum_{k=1}^N(Δt) J_k
    sum_J_k <- sum(rkou(N_Delta_t[i], p_up = p_up, eta_1 = eta_1, eta_2 = eta_2))

    # 4. Combine the drift, diffusion, and jump components to get the log-return
    # θ·Δt + σ·W(Δt) + \sum_{k=1}^N(Δt) J_k
    X[i] <- theta * Delta_t + diffusion + sum_J_k

    # 5. Compute pdf and cdf of ln(Xt/X0)
    f_X[i] <- dDoubleExpJumpDiff(X[i], Delta_t = Delta_t,
                                 mu = mu, sigma = sigma, lambda = lambda,
                                 p_up = p_up, eta_1 = eta_1, eta_2 = eta_2)

    F_X[i] <- pDoubleExpJumpDiff(X[i], Delta_t = Delta_t,
                                 mu = mu, sigma = sigma, lambda = lambda,
                                 p_up = p_up, eta_1 = eta_1, eta_2 = eta_2)

  }


  f_X_func <- (function(Delta_t, mu, sigma, lambda, p_up, eta_1, eta_2) {
    function(x) {
      sapply(x, function(xi) {
        dDoubleExpJumpDiff(xi, Delta_t = Delta_t,
                           mu = mu, sigma = sigma, lambda = lambda,
                           p_up = p_up, eta_1 = eta_1, eta_2 = eta_2)
      })
    }
  })(Delta_t, mu, sigma, lambda, p_up, eta_1, eta_2)

  F_X_func <- (function(Delta_t, mu, sigma, lambda, p_up, eta_1, eta_2) {
    function(x) {
      sapply(x, function(xi) {
        pDoubleExpJumpDiff(xi, Delta_t = Delta_t,
                           mu = mu, sigma = sigma, lambda = lambda,
                           p_up = p_up, eta_1 = eta_1, eta_2 = eta_2)
      })
    }
  })(Delta_t, mu, sigma, lambda, p_up, eta_1, eta_2)

  return(list(X = X,
              f_X = f_X, f_X_func = f_X_func,
              F_X = F_X, F_X_func = F_X_func))
}



# f_J(y) = p_{up} · η1e^{−η1y} 1[y ≥ 0] + (1-p_{up}) · η2e^{−η2y} 1[y < 0]
dKou <- function(x, p_up, eta_1, eta_2) {
  if (x < 0) {
    # For x < 0: use density eta_2 * exp(eta_2 * x)
    (1 - p_up) * eta_2 * exp(eta_2 * x)
  } else {
    # For x >= 0: use density eta_1 * exp(-eta_1 * x)
    p_up * eta_1 * exp(-eta_1 * x)
  }
}
# cdf
pKou <- function(x, p_up, eta_1, eta_2) {
  ifelse(x < 0,
         (1 - p_up) * exp(eta_2 * x),
         1 - p_up * exp(-eta_1 * x))
}
# quantile function
qKou <- function(u, p_up, eta_1, eta_2) {
  ifelse(u < (1 - p_up),
         (1/eta_2) * log(u / (1 - p_up)),
         - (1/eta_1) * log((1 - u) / p_up))
}
# sim
rkou <- function(n, p_up, eta_1, eta_2) {
  u <- runif(n)
  qKou(u, p_up = p_up, eta_1 = eta_1, eta_2 = eta_2)
}

# term1: (1 − λΔt)/(σ√Δt) \varphi ( (x - (μ - 0.5*σ^2 )Δt)/σ√Δt  )
# coef21: p_{up}η1 exp(0.5*η1^2σ^2Δt) * exp(-(x-(μ - 0.5*σ^2)Δt)η1)
# arg21: (x - (μ - 0.5*σ^2)Δt - η1*σ^2*Δt)/(σ√Δt)
# coef22: (1-p_{up})η2 exp(0.5*η2^2σ^2Δt) * exp((x-(μ - 0.5*σ^2)Δt)η2)
# arg22: (x - (μ - 0.5*σ^2)Δt - η2*σ^2*Δt)/(σ√Δt)
dDoubleExpJumpDiff <- function(x, Delta_t, mu, sigma, lambda, p_up, eta_1, eta_2) {

  coef1 <- (1 - lambda*Delta_t)/(sigma*sqrt(Delta_t))
  arg1 <- (x - (mu - 0.5*sigma^2)*Delta_t)/(sigma*sqrt(Delta_t))

  coef21 <- p_up*eta_1*exp(0.5*eta_1^2*sigma^2*Delta_t) * exp( -(x-(mu-0.5*sigma^2)*Delta_t)*eta_1 )
  coef22 <- (1-p_up)*eta_2 * exp(0.5*eta_2^2*sigma^2*Delta_t) * exp( (x-(mu-0.5*sigma^2)*Delta_t)*eta_2 )

  arg21 <- (x - (mu-0.5*sigma^2)*Delta_t - eta_1*sigma^2*Delta_t)/(sigma*sqrt(Delta_t))
  arg22 <- - (x - (mu-0.5*sigma^2)*Delta_t + eta_2*sigma^2*Delta_t)/(sigma*sqrt(Delta_t))


  return( coef1*dnorm(arg1) + lambda*Delta_t*(coef21*pnorm(arg21)+coef22*pnorm(arg22)) )

}

pDoubleExpJumpDiff <- function(x, Delta_t, mu, sigma, lambda, p_up, eta_1, eta_2) {

  b = 1/(sigma * sqrt(Delta_t))
  a_1 = b*(1 - lambda*Delta_t)

  a_21 = p_up*eta_1*exp(0.5*eta_1^2*sigma^2*Delta_t) * exp(eta_1*(mu-0.5*sigma^2)*Delta_t)
  a_22 = (1-p_up)*eta_2*exp(0.5*eta_2^2*sigma^2*Delta_t) * exp(-eta_2*(mu-0.5*sigma^2)*Delta_t)

  b = 1/(sigma * sqrt(Delta_t))

  c_1 = - b*(mu - 0.5*sigma^2)*Delta_t
  c_21 = c_1 - b*eta_1*sigma^2*Delta_t
  c_22 = c_1 + b*eta_2*sigma^2*Delta_t


  I_1 = (a_1/b)*pnorm(b*x +c_1)

  I_21 = -(a_21/eta_1) * ( exp(-eta_1*x) * pnorm(b*x+c_21)- exp(0.5*(eta_1/b)^2+eta_1*c_21/b) * pnorm(b*x+c_21+eta_1/b) )
  I_22 =  (a_22/eta_2) * ( exp(eta_2*x) * (1-pnorm(b*x+c_22)) + exp(0.5*(eta_2/b)^2-eta_2*c_22/b) * pnorm(b*x+c_22-eta_2/b) )

  return( I_1+lambda*Delta_t*(I_21+I_22) )

}


pDoubleExpJumpDiff_num <- function(x, Delta_t, mu, sigma, lambda, p_up, eta_1, eta_2) {

  coef1 <- (1 - lambda * Delta_t) / (sigma * sqrt(Delta_t))
  arg1  <- (x - (mu - 0.5 * sigma^2) * Delta_t) / (sigma * sqrt(Delta_t))

  # Integrand for I_1 with safety check
  integrand1 <- function(u) {
    coef21 <- p_up * eta_1 * exp(0.5 * eta_1^2 * sigma^2 * Delta_t) *
      exp(- (u - (mu - 0.5 * sigma^2) * Delta_t) * eta_1)
    arg21  <- (u - (mu - 0.5 * sigma^2) * Delta_t - eta_1 * sigma^2 * Delta_t) /
      (sigma * sqrt(Delta_t))
    res <- coef21 * pnorm(arg21)
    # Replace non-finite values with 0
    res[!is.finite(res)] <- 0
    return(res)
  }

  I_1 <- integrate(integrand1, lower = -Inf, upper = x)

  # Integrand for I_2 with safety check
  integrand2 <- function(u) {
    coef22 <- (1 - p_up) * eta_2 * exp(0.5 * eta_2^2 * sigma^2 * Delta_t) *
      exp((u - (mu - 0.5 * sigma^2) * Delta_t) * eta_2)
    arg22  <- - (u - (mu - 0.5 * sigma^2) * Delta_t + eta_2 * sigma^2 * Delta_t) /
      (sigma * sqrt(Delta_t))
    res <- coef22 * pnorm(arg22)
    res[!is.finite(res)] <- 0
    return(res)
  }

  I_2 <- integrate(integrand2, lower = -Inf, upper = x)

  # Extract the numerical values with I_1$value and I_2$value
  return(coef1 * pnorm(arg1) + lambda * Delta_t * (I_1$value + I_2$value))
}




# Delta_t = 1/4
# mu = 0.5
# sigma = 0.04
# lambda = 2
# p_up = 0.4
# eta_1 = 3
# eta_2 = 5
# for (x in seq(-3, 5, length.out = 1000)) {
#
#   val_num <- integrate(dDoubleExpJumpDiff,
#                        lower = -200,
#                        upper = x,
#                        Delta_t = Delta_t,
#                        mu = mu,
#                        sigma = sigma,
#                        lambda = lambda,
#                        p_up = p_up,
#                        eta_1 = eta_1,
#                        eta_2 = eta_2)$value
#
#
#
#
#   val <- pDoubleExpJumpDiff(x = x,
#                             Delta_t = Delta_t,
#                             mu = mu,
#                             sigma = sigma,
#                             lambda = lambda,
#                             p_up = p_up,
#                             eta_1 = eta_1,
#                             eta_2 = eta_2)
#
#   print(round(val_num-val,8))
# }









