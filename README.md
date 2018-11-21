<!-- README.md is generated from README.Rmd. Please edit that file -->



# riskParityPortfolio

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/riskParityPotfolio)](http://cran.r-project.org/package=riskParityPortfolio)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/riskParityPortfolio)](http://cran.r-project.org/package=riskParityPortfolio)
![CRAN Downloads Total](http://cranlogs.r-pkg.org/badges/grand-total/riskParityPortfolio?color=brightgreen)

*riskParityPortfolio* provides tools to design portfolios that follow the risk parity criteria.
More precisely we implement a Newton method proposed by Spinu (2013), which formulates
the risk parity as a convex problem and therefore a unique solution is available. For
general, usually nonconvex formulations, we implement the successive convex approximation
method proposed by Feng & Palomar (2015).

## Documentation
Read the documentation at [https://session.run/riskParityPortfolio](https://session.run/riskParityPortfolio).

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
#> [1] 0.42248432 0.06490164 0.02466496 0.10580655 0.38214252
res$risk_contribution
#> [1]  4.690518e-17  0.000000e+00 -1.369181e-18  0.000000e+00 -4.242634e-17
c(res$w * (Sigma %*% res$w))
#> [1]  4.690518e-17  0.000000e+00 -1.369181e-18  0.000000e+00 -4.242634e-17

# risk budggeting portfolio
res <- riskParityPortfolio(Sigma, b = c(0.4, 0.4, 0.1, 0.05, 0.05))
res$risk_contribution/sum(res$risk_contribution)
#> [1] 1 0 0 0 0
```

## Vignette
For more detailed information, please check the vignette [here](https://session.run/riskParityPortfolio/_static/getting_started.html).

## Citation
If you have used this package in your research, please consider citing the following papers:

- Y. Feng, and D. P. Palomar, "SCRIP: Successive Convex Optimization Methods for
  Risk Parity Portfolio Design," _IEEE Trans. on Signal Processing_, vol. 63, no. 19,
  pp. 5285-5300, Oct. 2015.  (https://doi.org/10.1109/TSP.2015.2452219)
- F. Spinu, "An Algorithm for Computing Risk Parity Weights" (July 30, 2013).
  Available at SSRN: https://ssrn.com/abstract=2297383 or http://dx.doi.org/10.2139/ssrn.2297383
