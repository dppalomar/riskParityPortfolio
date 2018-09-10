##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type="source")
# Installation from GitHub
devtools::install_github("dppalomar/riskParityPortfolio")
# Installation from CRAN
install.packages("spectralGraphTopology")
# Getting help
library(riskParityPortfolio)
help(package="riskParityPortfolio")
package?riskParityPortfolio
?riskParityPortfolioSCA


##
## Developer commands (http://r-pkgs.had.co.nz/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::install()
library(riskParityPortfolio)
#devtools::build()  # to generate the installation file

# Documentation
devtools::document()  #to generate all documentation via roxygen
?riskParityPortfolioSCA