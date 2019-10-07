## Changes in riskParityPortfolio version 0.2.1 (2019-10-07)

* A new section "A practical example using FAANG price data" was added to the vignette. This section is inspired by Tharsis Souza's blog post on risk parity: https://towardsdatascience.com/ray-dalio-etf-900edfe64b05


## Changes in riskParityPortfolio version 0.2.0 (2019-08-31)

* Included the R/Finance 2019 slides as an additional vignette.
* Included the slides on risk parity portfolio from the Convex Optimization course at 
  the Hong Kong Univ. of Science and Technology (HKUST) as an additional vignette.
* New plotting function implemented: barplotPortfolioRisk().
* General linear constraints now supported in the main function riskParityPortfolio()


## Changes in riskParityPortfolio version 0.1.2 (2019-06-01)

* Fixed some VignetteBuilder issues with CRAN.
* Refactored stopping criteria. [commit 350f622]
* Fixed bug where stocks names were being tossed out by C++ functions. [commit a02ffc4]


## Changes in riskParityPortfolio version 0.1.1 (2019-01-07)

* Revised vignette (fix name issue and include new section on algorithm description).
* Revise the error control of riskParityPortfolio().
* Check feasibility in riskParityPortfolio().
* Improved tests.


## Changes in riskParityPortfolio version 0.1.0 (2018-12-15)

* Initial release is on CRAN.
