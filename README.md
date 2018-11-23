<!-- README.md is generated from README.Rmd. Please edit that file -->



# riskParityPortfolio

[![codecov](https://codecov.io/gh/mirca/riskParityPortfolio/branch/master/graph/badge.svg)](https://codecov.io/gh/mirca/riskParityPortfolio)
[![Travis-CI-Badge](https://travis-ci.org/mirca/riskParityPortfolio.svg?branch=master)](https://travis-ci.org/mirca/riskParityPortfolio)
[![Build status](https://ci.appveyor.com/api/projects/status/dqjti1y461u7sjn8/branch/master?svg=true)](https://ci.appveyor.com/project/mirca/riskparityportfolio/branch/master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/riskParityPotfolio)](http://cran.r-project.org/package=riskParityPortfolio)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/riskParityPortfolio)](http://cran.r-project.org/package=riskParityPortfolio)
![CRAN Downloads Total](http://cranlogs.r-pkg.org/badges/grand-total/riskParityPortfolio?color=brightgreen)

*riskParityPortfolio* provides tools to design portfolios that follow the risk parity criteria.
More precisely we implement a Newton method proposed by Spinu (2013), which formulates
the risk parity as a convex problem and therefore a unique solution is available. For
general, usually nonconvex formulations, we implement the successive convex approximation
method proposed by Feng & Palomar (2015).


## Installation

```r
# Installation from GitHub
install.packages("devtools")
devtools::install_github("dppalomar/riskParityPortfolio")

# Getting help
library(riskParityPortfolio)
help(package = "riskParityPortfolio")
package?riskParityPortfolio
?riskParityPortfolio

# Citing this work
citation("riskParityPortfolio")
```

## Usage of `riskParityPortfolio`

```r
library(riskParityPortfolio)

# create covariance matrix
N <- 5
V <- matrix(rnorm(N^2), nrow = N)
Sigma <- cov(V)

# risk-parity portfolio
res <- riskParityPortfolio(Sigma)
names(res)
#> [1] "w"                 "risk_contribution"
res$w
#> [1] 0.12642816 0.02512529 0.02856941 0.54993677 0.26994037
res$risk_contribution
#> [1] 0.0009458898 0.0009458898 0.0009458898 0.0009458898 0.0009458898
c(res$w * (Sigma %*% res$w))
#> [1] 0.0009458898 0.0009458898 0.0009458898 0.0009458898 0.0009458898

# risk budggeting portfolio
res <- riskParityPortfolio(Sigma, b = c(0.4, 0.4, 0.1, 0.05, 0.05))
res$risk_contribution/sum(res$risk_contribution)
#> [1] 0.40 0.40 0.10 0.05 0.05
```

## Documentation
For more detailed information, please check the vignette
[here](https://htmlpreview.github.io/?https://github.com/dppalomar/riskParityPortfolio/blob/master/vignettes/RiskParityPortfolio-vignette.html)
or the package webpage [https://mirca.github.io/riskParityPortfolio](https://mirca.github.io/riskParityPortfolio).

## Citation
If you have used this package in your research, please consider citing the following papers:

- Y. Feng, and D. P. Palomar, "SCRIP: Successive Convex Optimization Methods for
  Risk Parity Portfolio Design," _IEEE Trans. on Signal Processing_, vol. 63, no. 19,
  pp. 5285-5300, Oct. 2015.  (https://doi.org/10.1109/TSP.2015.2452219)
- F. Spinu, "An Algorithm for Computing Risk Parity Weights" (July 30, 2013).
  Available at SSRN: https://ssrn.com/abstract=2297383 or http://dx.doi.org/10.2139/ssrn.2297383
  
## Links
Package: [GitHub](https://github.com/dppalomar/riskParityPortfolio).
README file: [GitHub-readme](https://htmlpreview.github.io/?https://github.com/dppalomar/riskParityPortfolio/blob/master/README.html).
Vignette: [GitHub-html-vignette](https://htmlpreview.github.io/?https://github.com/dppalomar/riskParityPortfolio/blob/master/vignettes/RiskParityPortfolio-vignette.html) and
[GitHub-pdf-vignette](https://docs.google.com/viewer?url=https://github.com/dppalomar/riskParityPortfolio/raw/master/vignettes/RiskParityPortfolio-vignette.pdf)

