<!-- README.md is generated from README.Rmd. Please edit that file -->
riskParityPortfolio
===================

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/riskParityPotfolio)](http://cran.r-project.org/package=riskParityPortfolio)
[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/riskParityPortfolio)](http://cran.r-project.org/package=riskParityPortfolio)
![CRAN Downloads
Total](http://cranlogs.r-pkg.org/badges/grand-total/riskParityPortfolio?color=brightgreen)

This package provides tools to design portfolios that follow the risk
parity criteria. More precisely we implement a Newton method proposed by
Spinu (2013), which formulates the risk parity as a convex problem and
therefore a unique solution is available. For general, usually nonconvex
formulations, we implement the successive convex approximation method
proposed by Feng & Palomar (2015).

This package is based on the papers:

- Y. Feng, and D. P. Palomar,
“SCRIP: Successive Convex Optimization Methods for Risk Parity Portfolio
Design” *IEEE Trans. on Signal Processing*, vol. 63, no. 19,
pp. 5285-5300, Oct. 2015. <doi:10.1109/TSP.2015.2452219>
- F. Spinu, “An Algorithm for Computing Risk Parity Weights” (July 30, 2013). Available
at SSRN: <https://ssrn.com/abstract=2297383> or
<http://dx.doi.org/10.2139/ssrn.2297383>

Installation
------------

``` r
# Installation from GitHub
install.packages("devtools")
devtools::install_github("dppalomar/riskParityPortfolio")

# Getting help
library(riskParityPortfolio)
help(package = "riskParityPortfolio")
package?riskParityPortfolio
?spIndexTrack

# Citing this work
citation("riskParityPortfolio")
```

Vignette
--------

For more detailed information, please check the vignette:
[GitHub-html-vignette](https://rawgit.com/dppalomar/riskParityPortfolio/master/vignettes/RiskParityPortfolio-vignette.html),
[GitHub-pdf-vignette](https://rawgit.com/dppalomar/riskParityPortfolio/master/vignettes/RiskParityPortfolio-vignette.pdf),
[CRAN-pdf-vignette](https://cran.r-project.org/web/packages/riskParityPortfolio/vignettes/RiskParityPortfolio-vignette.pdf).

Usage of `riskParityPortfolio`
------------------------------

``` r
library(riskParityPortfolio)

# Long-only portfolio design, Sigma is the covariance matrix.
portfolio_long <- riskParityPortfolioNewton(Sigma)
# Long/short portfolio design using the successive convex approximation method.
portfolio_sca <- riskParityPortfolioSCA(Sigma)
# Long/short portfolio design using the general non-linear solvers.
portfolio_gen <- riskParityPortfolioGenSolver(Sigma)
```

Citation
--------

If you have used this package in your research, please consider citing
the following papers:

-   Y. Feng, and D. P. Palomar, “SCRIP: Successive Convex Optimization
    Methods for Risk Parity Portfolio Design” *IEEE Trans. on Signal
    Processing*, vol. 63, no. 19, pp. 5285-5300, Oct. 2015.
    <doi:10.1109/TSP.2015.2452219>
-   F. Spinu, “An Algorithm for Computing Risk Parity Weights” (July 30,
    2013). Available at SSRN: <https://ssrn.com/abstract=2297383> or
    <http://dx.doi.org/10.2139/ssrn.2297383>

Links
-----

Package: [CRAN](https://cran.r-project.org/package=riskParityPortfolio)
and [GitHub](https://github.com/dppalomar/riskParityPortfolio).

README file:
[GitHub-readme](https://rawgit.com/dppalomar/riskParityPortfolio/master/README.html)
and
[CRAN-readme](https://cran.r-project.org/web/packages/riskParityPortfolio/README.html).

Vignette:
[GitHub-html-vignette](https://rawgit.com/dppalomar/riskParityPortfolio/master/vignettes/RiskParityPortfolio-vignette.html)
and
[GitHub-pdf-vignette](https://rawgit.com/dppalomar/riskParityPortfolio/master/vignettes/RiskParityPortfolio-vignette.pdf),
[CRAN-pdf-vignette](https://cran.r-project.org/web/packages/riskParityPortfolio/vignettes/RiskParityPortfolio-vignette.pdf).
