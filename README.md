# riskParityPortfolio

[![codecov](https://codecov.io/gh/mirca/riskParityPortfolio/branch/master/graph/badge.svg?token=NhOMM6MOAv)](https://codecov.io/gh/mirca/riskParityPortfolio)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/riskParityPortfolio)](https://CRAN.R-project.org/package=riskParityPortfolio)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/riskParityPortfolio)](https://CRAN.R-project.org/package=riskParityPortfolio)
[![CRAN Downloads Total](https://cranlogs.r-pkg.org/badges/grand-total/riskParityPortfolio?color=brightgreen)](https://CRAN.R-project.org/package=riskParityPortfolio)
[![Rcpp](https://img.shields.io/badge/powered%20by-Rcpp-orange.svg?style=flat)](http://www.rcpp.org/)


**riskParityPortfolio** provides tools to design risk parity portfolios.
In its simplest form, we consider the convex formulation with a unique solution proposed by
[Spinu (2013)](https://dx.doi.org/10.2139/ssrn.2297383) and use a cyclical method inspired by
[Griveau-Billion (2013)](https://arxiv.org/pdf/1311.4057.pdf). For more general formulations,
which are usually nonconvex, we implement the successive convex approximation
method proposed by [Feng & Palomar (2015)](https://doi.org/10.1109/TSP.2015.2452219).

The latest stable version of **riskParityPortfolio** is available at https://CRAN.R-project.org/package=riskParityPortfolio.

The latest development version of **riskParityPortfolio** is available at https://github.com/dppalomar/riskParityPortfolio.

**Check out the documentation here: [https://mirca.github.io/riskParityPortfolio](https://mirca.github.io/riskParityPortfolio/).**

## Installation
To install the latest stable version of **riskParityPortfolio** from CRAN, run the following commands in R:

```r
> install.packages("riskParityPortfolio")
```

To install the development version of **riskParityPortfolio** from GitHub, run the following commands in R:

```r
> install.packages("devtools")
> devtools::install_github("dppalomar/riskParityPortfolio")
```

To get help:

```r
> library(riskParityPortfolio)
> help(package = "riskParityPortfolio")
> package?riskParityPortfolio
> ?riskParityPortfolio
```

Please cite **riskParityPortfolio** in publications:

```r
> citation("riskParityPortfolio")
```

You can also get **riskParityPortfolio** from Docker as follows:
```
$ docker pull mirca/riskparityportfolio
```

#### Microsoft Windows
On MS Windows environments, make sure to install the most recent version of
**Rtools**.

### Python

A Python3 implementation of this package is currently under development at [https://github.com/dppalomar/riskparity.py](https://github.com/dppalomar/riskparity.py).
Its stable version is available in PYPI and can be installed as follows:
```
$ pip install riskparityportfolio
```

Alternatively, the development version can be installed as
```
$ git clone https://github.com/dppalomar/riskparity.py
$ cd riskparity.py
$ pip install -e .
```

### Rust

A Rust version is avaiable at: [https://github.com/mirca/riskparity.rs](https://github.com/mirca/riskparity.rs).

## Basic usage


```r
library(riskParityPortfolio)

set.seed(42)
# create covariance matrix
N <- 5
V <- matrix(rnorm(N^2), ncol = N)
Sigma <- cov(V)

# risk parity portfolio
res <- riskParityPortfolio(Sigma)
names(res)
#> [1] "w"                          "relative_risk_contribution"
#> [3] "obj_fun"                    "is_feasible"
res$w
#> [1] 0.32715962 0.27110678 0.14480081 0.09766356 0.15926922

# risk budgeting portfolio
res <- riskParityPortfolio(Sigma, b = c(0.4, 0.4, 0.1, 0.05, 0.05))
res$relative_risk_contribution
#> [1] 0.40 0.40 0.10 0.05 0.05
```

## Documentation
For more detailed information, please check the
[vignette](https://CRAN.R-project.org/package=riskParityPortfolio/vignettes/RiskParityPortfolio.html).

## Citation
If you find this package useful in your research, please consider citing the following works:

- J. V. de M. Cardoso and D. P. Palomar (2019). riskParityPortfolio:
  Design of Risk Parity Portfolios. R package version 0.2.1.
  <https://CRAN.R-project.org/package=riskParityPortfolio>
- Y. Feng, and D. P. Palomar (2015). SCRIP: Successive Convex Optimization Methods for
  Risk Parity Portfolio Design. _IEEE Trans. on Signal Processing_, vol. 63, no. 19,
  pp. 5285-5300. <https://doi.org/10.1109/TSP.2015.2452219>
- F. Spinu (2013). An Algorithm for Computing Risk Parity Weights.
  <https://dx.doi.org/10.2139/ssrn.2297383>
- T. Griveau-Billion, J. Richard, and T. Roncalli (2013). A fast algorithm for computing High-dimensional risk parity portfolios. <https://arxiv.org/pdf/1311.4057.pdf>


## Contributing
We welcome all sorts of contributions. Please feel free to open an issue
to report a bug or discuss a feature request.


## Links
Package: [CRAN](https://CRAN.R-project.org/package=riskParityPortfolio) and [GitHub](https://github.com/dppalomar/riskParityPortfolio).

README file: [GitHub-readme](https://github.com/dppalomar/riskParityPortfolio/blob/master/README.md).

Vignettes: [CRAN-vignette](https://CRAN.R-project.org/package=riskParityPortfolio/vignettes/RiskParityPortfolio.html),
[slides R/Finance 2019](https://docs.google.com/viewer?url=https://github.com/dppalomar/riskParityPortfolio/raw/master/vignettes/RFinance2019-slides.pdf),
[slides RPP - Convex Optimization Course (HKUST)](https://docs.google.com/viewer?url=https://github.com/dppalomar/riskParityPortfolio/raw/master/vignettes/slides-ConvexOptimizationCourseHKUST.pdf),
[slides HKML meetup 2020](https://speakerdeck.com/mirca/breaking-down-risk-parity-portfolios-a-practical-open-source-implementation), and
[tutorial - Data-driven Portfolio Optimization Course (HKUST)](https://www.youtube.com/watch?v=xb1Xxf5LQks)


## Disclaimer
The information, software, and any additional resources contained in this repository are not intended as,
and shall not be understood or construed as, financial advice.
Past performance is not a reliable indicator of future results and investors may not recover the full
amount invested.
The [authors](https://github.com/dppalomar/riskParityPortfolio/blob/master/AUTHORS.md) of this repository
accept no liability whatsoever for any loss or damage you may incur.  Any opinions expressed in this repository
are from the personal research and experience of the [authors](https://github.com/dppalomar/riskParityPortfolio/blob/master/AUTHORS.md) and are intended as educational material.

