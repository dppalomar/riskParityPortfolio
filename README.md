<!-- README.md is generated from README.Rmd. Please edit that file -->
riskParityPortfolio
===================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/riskParityPortfolio)](https://cran.r-project.org/package=riskParityPortfolio)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/riskParityPortfolio)](https://cran.r-project.org/package=riskParityPortfolio)
![CRAN Downloads
Total](https://cranlogs.r-pkg.org/badges/grand-total/riskParityPortfolio?color=brightgreen)

[![PyPI
version](https://badge.fury.io/py/riskparityportfolio.svg)](https://badge.fury.io/py/riskparityportfolio)
[![Downloads](https://pepy.tech/badge/riskparityportfolio)](https://pepy.tech/project/riskparityportfolio)

[![codecov](https://codecov.io/gh/mirca/riskParityPortfolio/branch/master/graph/badge.svg)](https://codecov.io/gh/mirca/riskParityPortfolio)
[![Travis-CI-Badge](https://travis-ci.org/mirca/riskParityPortfolio.svg?branch=master)](https://travis-ci.org/mirca/riskParityPortfolio)
[![Build
status](https://ci.appveyor.com/api/projects/status/dqjti1y461u7sjn8/branch/master?svg=true)](https://ci.appveyor.com/project/mirca/riskparityportfolio/branch/master)
[![CircleCI](https://circleci.com/gh/mirca/riskParityPortfolio.svg?style=svg)](https://circleci.com/gh/mirca/riskParityPortfolio)
[![Docker Build
Status](https://img.shields.io/docker/build/mirca/riskparityportfolio.svg)](https://hub.docker.com/r/mirca/riskparityportfolio/)
[![Build
Status](https://dev.azure.com/jvmirca/riskParityPortfolio/_apis/build/status/mirca.riskParityPortfolio?branchName=master)](https://dev.azure.com/jvmirca/riskParityPortfolio/_build/latest?definitionId=1&branchName=master)

**riskParityPortfolio** provides tools to design risk parity portfolios.
In its simplest form, we consider the convex formulation with a unique
solution proposed by [Spinu
(2013)](https://dx.doi.org/10.2139/ssrn.2297383) and use a cyclical
method inspired by [Griveau-Billion
(2013)](https://arxiv.org/pdf/1311.4057.pdf). For more general
formulations, which are usually nonconvex, we implement the successive
convex approximation method proposed by [Feng & Palomar
(2015)](https://doi.org/10.1109/TSP.2015.2452219).

The latest stable version of **riskParityPortfolio** is available at
<https://CRAN.R-project.org/package=riskParityPortfolio>.

The latest development version of **riskParityPortfolio** is available
at <https://github.com/dppalomar/riskParityPortfolio>.

**Check out the
[documentation](https://mirca.github.io/riskParityPortfolio).**

Installation
------------

To install the latest stable version of **riskParityPortfolio**, run the
following commands in R:

``` r
install.packages("riskParityPortfolio")
```

To install the development version of **riskParityPortfolio**, run the
following commands in R:

``` r
install.packages("devtools")
devtools::install_github("dppalomar/riskParityPortfolio")
```

To get help:

``` r
library(riskParityPortfolio)
help(package = "riskParityPortfolio")
package?riskParityPortfolio
?riskParityPortfolio
```

Please cite **riskParityPortfolio** in publications:

``` r
citation("riskParityPortfolio")
```

You can also get **riskParityPortfolio** from Docker as follows:

    docker pull mirca/riskparityportfolio

#### Microsoft Windows

On MS Windows environments, make sure to install the most recent version
of **Rtools**.

### Python

A Python3 implementation of the vanilla method is available in PYPI and
can be installed as follows:

    pip install riskparityportfolio

Usage of **riskParityPortfolio**
--------------------------------

### R

``` r
library(riskParityPortfolio)

set.seed(42)
# create covariance matrix
N <- 5
V <- matrix(rnorm(N^2), ncol = N)
Sigma <- cov(V)

# risk parity portfolio
res <- riskParityPortfolio(Sigma)
names(res)
#> [1] "w"                 "risk_contribution"
res$w
#> [1] 0.32715962 0.27110678 0.14480081 0.09766356 0.15926922
res$risk_contribution
#> [1] 0.03857039 0.03857039 0.03857039 0.03857039 0.03857039
c(res$w * (Sigma %*% res$w))
#> [1] 0.03857039 0.03857039 0.03857039 0.03857039 0.03857039

# risk budggeting portfolio
res <- riskParityPortfolio(Sigma, b = c(0.4, 0.4, 0.1, 0.05, 0.05))
res$risk_contribution/sum(res$risk_contribution)
#> [1] 0.40 0.40 0.10 0.05 0.05
```

### Python

``` python
import numpy as np
import riskparityportfolio as rpp
np.random.seed(42)
# creates a correlation matrix from time-series of five assets
x = np.random.normal(size=1000).reshape((5, -1))
corr = x @ x.T
# create the desired risk budgeting vector
b = np.ones(len(corr)) / len(corr)
# design the portfolio
w = rpp.design(corr, b)
print(w)
# compute the risk budgeting
#> [0.21075375 0.21402865 0.20205399 0.16994639 0.20321721]
rc = w @ (corr * w)
print(rc / np.sum(rc))
# let's try a different budget
#> [0.2 0.2 0.2 0.2 0.2]
b = np.array([0.01, 0.09, .1, .1, .7])
w = rpp.design(corr, b)
print(w)
#> [0.06178354 0.19655744 0.16217134 0.12808275 0.45140493]
rc = w @ (corr * w)
print(rc / np.sum(rc))
#> [0.01 0.09 0.1  0.1  0.7 ]
```

Documentation
-------------

For more detailed information, please check the [CRAN
vignette](https://CRAN.R-project.org/package=riskParityPortfolio/vignettes/RiskParityPortfolio.html)
and [GitHub
vignette](https://raw.githack.com/dppalomar/riskParityPortfolio/master/vignettes/RiskParityPortfolio.html).

Citation
--------

If you find this package useful in your research, please consider citing
the following works:

-   J. V. de M. Cardoso and D. P. Palomar (2019). riskParityPortfolio:
    Design of Risk Parity Portfolios. R package version 0.1.2.
    <https://CRAN.R-project.org/package=riskParityPortfolio>
-   Y. Feng, and D. P. Palomar (2015). SCRIP: Successive Convex
    Optimization Methods for Risk Parity Portfolio Design. *IEEE Trans.
    on Signal Processing*, vol. 63, no. 19, pp. 5285-5300.
    <https://doi.org/10.1109/TSP.2015.2452219>
-   F. Spinu (2013). An Algorithm for Computing Risk Parity Weights.
    <https://dx.doi.org/10.2139/ssrn.2297383>
-   T. Griveau-Billion, J. Richard, and T. Roncalli (2013). A fast
    algorithm for computing High-dimensional risk parity portfolios.
    <https://arxiv.org/pdf/1311.4057.pdf>

Contributing
------------

We welcome all sorts of contributions. Please feel free to open an issue
to report a bug or discuss a feature request.

Links
-----

Package: [CRAN](https://CRAN.R-project.org/package=riskParityPortfolio)
and [GitHub](https://github.com/dppalomar/riskParityPortfolio).

README file:
[CRAN-readme](https://CRAN.R-project.org/package=riskParityPortfolio/readme/README.html)
and
[GitHub-readme](https://github.com/dppalomar/riskParityPortfolio/blob/master/README.md).

Vignette:
[CRAN-html-vignette](https://CRAN.R-project.org/package=riskParityPortfolio/vignettes/RiskParityPortfolio.html),
[CRAN-pdf-vignette](https://CRAN.R-project.org/package=riskParityPortfolio/vignettes/RiskParityPortfolio-pdf.pdf),
[GitHub-html-vignette](https://raw.githack.com/dppalomar/riskParityPortfolio/master/vignettes/RiskParityPortfolio.html),
and
[GitHub-pdf-vignette](https://docs.google.com/viewer?url=https://github.com/dppalomar/riskParityPortfolio/raw/master/vignettes/RiskParityPortfolio-pdf.pdf).
