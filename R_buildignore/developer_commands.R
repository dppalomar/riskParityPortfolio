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
citation("riskParityPortfolio")
vignette(package = "riskParityPortfolio")


##
## Developer commands (http://r-pkgs.had.co.nz/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::install()
#devtools::install(build_vignettes = TRUE)
library(riskParityPortfolio)

# Documentation
devtools::document()  # to generate all documentation via roxygen
?riskParityPortfolio


# Code tests
#devtools::use_testthat()  # the first time
devtools::test()
#covr::package_coverage()  # coverage of tests
#goodpractice::gp()  # overall checks


# CRAN check and submission (http://r-pkgs.had.co.nz/release.html)
#  checklist: https://kalimu.github.io/post/checklist-for-r-package-submission-to-cran/
devtools::check()
rcmdcheck::rcmdcheck()
devtools::build()
#devtools::build_win()  #to check under windows
#R CMD build .  # this is to generate tarball
#R CMD build . --compact-vignettes=gs+qpdf  
#R CMD check riskParityPortfolio_0.1.0.9000.tar.gz --as-cran  # this is before submission to CRAN
#R CMD install riskParityPortfolio_0.1.0.9000.tar.gz
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html
