# ddfs

The **ddfs** R package provides nonparametric estimation of univariate and
bivariate **pdfs** and **cdfs** using variable-knot splines. This method leverages
an adaptive data-driven knot selection process (from Geometrically Designed
Splines) and enforces necessary shape constraints (non-negativity and monotonicity)
through constrained spline estimation, combining theoretical robustness with
practical usability.
---

## Installation

```r
devtools::install_github("emilioluissaenzguillen/ddfs")
```
```r
## Example
library(ddfs)

# Simulate data and fit a CDF
set.seed(123)
sim <- sim.dist(N = 1000, ex = "Gaussian")
ddfs_fit <- ddfs(sim$X)

# Plot estimated pdf/cdf
plot(ddfs_fit, type = "density", f = sim$f_X_func)
plot(ddfs_fit, type = "distribution", f = sim$F_X_func)
```

License:
========

This package is free and open source software, licensed under GPL-3

