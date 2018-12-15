##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type="source")
# Installation from GitHub
devtools::install_github("dppalomar/riskParityPortfolio")
# Installation from CRAN
install.packages("riskParityPortfolio")
# Getting help
library(riskParityPortfolio)
help(package = "riskParityPortfolio")
package?riskParityPortfolio
?riskParityPortfolio


##
## Developer commands (http://r-pkgs.had.co.nz/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::install()
library(riskParityPortfolio)
#devtools::build()  # to generate the installation file

# Documentation
devtools::document()  #to generate all documentation via roxygen
?riskParityPortfolio


# Code tests
#devtools::use_testthat()  # the first time
devtools::test()
#covr::package_coverage()  #coverage of tests
#goodpractice::gp()  # overall checks


# CRAN check and submission
devtools::check()
rcmdcheck::rcmdcheck()
devtools::build()
#R CMD build .
#R CMD build . --compact-vignettes=gs+qpdf  # this is to generate tarball
#R CMD check riskParityPortfolio_0.1.0.tar.gz --as-cran  # this is before submission to CRAN
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html

# An alternative is to upload to CRAN via devtools:
#devtools::build_win()  #to check under windows
#devtools::release(args = "--compact-vignettes=gs+qpdf")  #for CRAN

