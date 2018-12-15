---
output:
  html_document:
    variant: markdown_github
    keep_md: true
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# riskParityPortfolio




The package *riskParityPortfolio* provides tools to design risk-parity 
portfolios. In its simplest form, we consider the convex formulation 
with a unique solution proposed by Spinu (2013) and use a cyclical 
method inspired by Griveau-Billion (2013). For more general formulations, 
which are usually nonconvex, we implement the successive convex approximation 
method proposed by Feng & Palomar (2015).


## Installation
To install ``riskParityPortfolio`` type the following inside an R session:


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

You can also get ``riskParityPortfolio`` from Docker as follows:
```
docker pull mirca/riskparityportfolio
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
#> [1] 0.08255655 0.47600209 0.16441764 0.10521082 0.17181289
res$risk_contribution
#> [1] 0.03587918 0.03587918 0.03587918 0.03587918 0.03587918
c(res$w * (Sigma %*% res$w))
#> [1] 0.03587918 0.03587918 0.03587918 0.03587918 0.03587918

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
If you have used this package in your research, please consider citing the following sources:

- J. V. de M. Cardoso and D. P. Palomar (2018). riskParityPortfolio:
  Design of Risk Parity Portfolios. R package version 0.1.0.
  https://CRAN.R-project.org/package=riskParityPortfolio
- Y. Feng, and D. P. Palomar, "SCRIP: Successive Convex Optimization Methods for
  Risk Parity Portfolio Design," _IEEE Trans. on Signal Processing_, vol. 63, no. 19,
  pp. 5285-5300, Oct. 2015.  (https://doi.org/10.1109/TSP.2015.2452219)
- F. Spinu, "An Algorithm for Computing Risk Parity Weights," 2013.
  Available at SSRN: https://ssrn.com/abstract=2297383 or http://dx.doi.org/10.2139/ssrn.2297383
- T. Griveau-Billion, J. Richard, and T. Roncalli, "A fast algorithm for computing High-dimensional risk parity portfolios," 2013.
  ArXiv preprint: https://arxiv.org/pdf/1311.4057.pdf

## Links
Package: [GitHub](https://github.com/dppalomar/riskParityPortfolio).
README file: [GitHub-readme](https://htmlpreview.github.io/?https://github.com/dppalomar/riskParityPortfolio/blob/master/README.html).
Vignette: [GitHub-html-vignette](https://htmlpreview.github.io/?https://github.com/dppalomar/riskParityPortfolio/blob/master/vignettes/RiskParityPortfolio-vignette.html) and
[GitHub-pdf-vignette](https://docs.google.com/viewer?url=https://github.com/dppalomar/riskParityPortfolio/raw/master/vignettes/RiskParityPortfolio-vignette.pdf)

