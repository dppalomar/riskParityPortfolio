<!-- README.md is generated from README.Rmd. Please edit that file -->



# riskParityPortfolio

This package provides tools to design portfolios that follow the risk parity criteria.
More precisely we implement a Newton method proposed by Spinu (2013), which formulates
the risk parity as a convex problem and therefore a unique solution is available. For
general, usually nonconvex formulations, we implement the successive convex approximation
method proposed by Feng & Palomar (2015).

This package is based on the papers:
- Y. Feng, and D. P. Palomar, "SCRIP: Successive Convex Optimization Methods for
  Risk Parity Portfolio Design" _IEEE Trans. on Signal Processing_, vol. 63, no. 19,
  pp. 5285-5300, Oct. 2015. (https://doi.org/10.1109/TSP.2015.2452219)
- F. Spinu, "An Algorithm for Computing Risk Parity Weights" (July 30, 2013).
  Available at SSRN: https://ssrn.com/abstract=2297383 or http://dx.doi.org/10.2139/ssrn.2297383

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


## Vignette
For more detailed information, please check the vignette:
[GitHub-html-vignette](https://rawgit.com/dppalomar/riskParityPortfolio/master/vignettes/RiskParityPortfolio-vignette.html),
[GitHub-pdf-vignette](https://rawgit.com/dppalomar/riskParityPortfolio/master/vignettes/RiskParityPortfolio-vignette.pdf).


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
res$w
res$risk_contribution
c(res$w * (Sigma %*% res$w))

# risk budggeting portfolio
res <- riskParityPortfolio(Sigma, b = c(0.4, 0.4, 0.1, 0.05, 0.05))
res$risk_contribution/sum(res$risk_contribution)
```

## Citation
If you have used this package in your research, please consider citing the following papers:

- Y. Feng, and D. P. Palomar, "SCRIP: Successive Convex Optimization Methods for
  Risk Parity Portfolio Design" _IEEE Trans. on Signal Processing_, vol. 63, no. 19,
  pp. 5285-5300, Oct. 2015.  (https://doi.org/10.1109/TSP.2015.2452219)
- F. Spinu, "An Algorithm for Computing Risk Parity Weights" (July 30, 2013).
  Available at SSRN: https://ssrn.com/abstract=2297383 or http://dx.doi.org/10.2139/ssrn.2297383
