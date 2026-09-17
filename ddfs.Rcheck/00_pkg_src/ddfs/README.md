# ddfs

The **ddfs** R package provides nonparametric estimation of univariate and
multivariate **pdfs** and **cdfs** using variable-knot splines.

This method leverages an adaptive data-driven knot selection process (borrowed 
from Geometrically Designed Splines) and enforces necessary constraints
(non-negativity and monotonicity) through constrained spline estimation, combining
theoretical robustness with practical usability.

---

## Installation

```r
devtools::install_github("emilioluissaenzguillen/ddfs")
```

## Example

```r
library(ddfs)

# Simulate data and fit a CDF
set.seed(123)
sim <- sim.dist(N = 1000, ex = "Gaussian")
ddfs_fit <- ddfs(sim$X)

# Plot estimated pdf/cdf
plot(ddfs_fit, fit = "pdf", f = sim$f_X_func)
plot(ddfs_fit, fit = "cdf", f = sim$F_X_func)
```

## Higher-dimensional fits (development)

Fits with two or more dimensions use the multidimensional GeDS engine and
require a GeDS build that includes it. Above two dimensions, missing tuning
parameters use provisional fixed defaults: `min.intknots = 10`, `beta = 0.1`,
`phi_F = 0.95`, and `q_F = 2`, with density order `n = 4`. Explicit values override
these defaults. They are based on preliminary 3D Gaussian experiments, not yet
validated across higher-dimensional distributions. Automatic selection remains
available for one or two dimensions. Order 4 can be costly in higher dimensions.

```r
set.seed(123)
XYZ <- cbind(x = runif(100), y = runif(100), z = runif(100))
fit <- ddfs(XYZ, n = 3, min.intknots = 1, max.intknots = 3,
            beta = 0.1, phi_F = 0.95, q_F = 1)
d.ddfs(fit, XYZ)
p.ddfs(fit, XYZ)
knots(fit, options = "internal")
```

Named prediction columns are matched to the training names; unnamed columns must
use the training order. These fits support PDF/CDF prediction, coefficient and
knot extraction, and summaries. Direct plotting, quantiles, and random sampling
are not supported. Only `stoptype = "RDMD"` and `schoenberg = FALSE` are supported.
The `max.coef` limit guards tensor-basis size, which grows rapidly with dimension.
Univariate fitting is unchanged. Bivariate fitting now uses the same internal
engine as higher-dimensional fitting, while retaining its existing public object
layout and methods. The legacy bivariate fitter remains a regression-test reference.

### Computational limits

With density order `n` and `k[j]` internal knots on coordinate `j`, the PDF
has `prod(n + k)` coefficients and the CDF has `prod(n + 1 + k)`.
The dense training matrices alone therefore require approximately
`8 * N * (prod(n + k) + prod(n + 1 + k))` bytes. This excludes temporary
copies, knot selection, and optimization; it is **not** peak memory usage.
`max.coef` limits coefficient counts, not total memory, and prediction memory
also grows with the number of query rows. The empirical CDF currently costs
quadratic time in the number of observations (times the dimension).

A small local benchmark (2026-09-03, Windows, R 4.6.1, development GeDS 0.3.5)
gave the following single-run results for cubic densities (`n = 4`), 100
observations, and three fitting iterations:

| Dimensions | Elapsed seconds | PDF coefficients | CDF coefficients | Training matrices, MiB |
| --- | ---: | ---: | ---: | ---: |
| 3 | 0.64 | 100 | 180 | 0.21 |
| 4 | 2.20 | 400 | 900 | 0.99 |
| 5 | 18.15 | 1,600 | 4,500 | 4.65 |

These are illustrative timings, not performance guarantees. For each dimension
`d`, the benchmark used `set.seed(910 + d)`, `matrix(runif(100 * d), ncol = d)`,
then scaled coordinate `j` by `j` and shifted it by `-j`. Fitting used
`min.intknots = max.intknots = 2`, `q_F = 1`, `beta = 0.1`, `phi_F = 0.95`,
and `tails_count_threshold = 0`; package loading was excluded from timing.

The 4D/5D regression tests cover density orders 2–4, normalization by independent
quadrature, joint-CDF rectangle probabilities, predictions, and coefficient limits.
The initial benchmark exposed integer-overflow warnings in GeDS's
`detectTensorMeshND()` for the 5D scattered samples. The local GeDS source now
checks whether a complete grid is possible before constructing cell identifiers,
avoiding that overflow while preserving the fast solver for complete grids.
Use the updated GeDS source to obtain this fix; the benchmark timings above
precede the fix.

License:
========

This package is free and open source software, licensed under GPL-3

