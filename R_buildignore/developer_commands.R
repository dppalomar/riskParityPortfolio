##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type = "source")
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
#tools::showNonASCIIfile("R/riskParityPortfolio.R")


##
## Developer commands (https://r-pkgs.org/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::document()  #or Ctrl-Shift-L (to generate all documentation via roxygen)
devtools::install(build_vignettes = TRUE)
library(riskParityPortfolio)
#tools::showNonASCIIfile("fit_mvt.R")


# Code tests (https://codecov.io/gh/mirca/riskParityPortfolio)
devtools::test()
#covr::package_coverage()  # coverage of tests


# Compress pdf vignettes
#tools::compactPDF("vignettes/slides-RFinance2019.pdf", gs_quality = "ebook")
#tools::compactPDF("vignettes/slides-ConvexOptimizationCourseHKUST.pdf", gs_quality = "ebook")


# CRAN check and submission (https://r-pkgs.org/release.html)
#  checklist: https://kalimu.github.io/post/checklist-for-r-package-submission-to-cran/
devtools::check()  # run_dont_test = TRUE
rcmdcheck::rcmdcheck()  # build_args = "--run-donttest"
devtools::build()
#devtools::check_win_release()  #to check under windows
#R CMD build .  # this is to generate tarball
#R CMD build . --compact-vignettes=gs+qpdf  
#R CMD check riskParityPortfolio_0.2.2.tar.gz --as-cran --run-donttest  # this is before submission to CRAN
#R CMD install riskParityPortfolio_0.2.2.tar.gz
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html



## Reverse dependencies
# tools::dependsOnPkgs("riskParityPortfolio")
# tools::package_dependencies("riskParityPortfolio", reverse = TRUE)
# devtools::revdep(pkg = "riskParityPortfolio") 
# revdepcheck::revdep_check()  #devtools::install_github("r-lib/revdepcheck")
# remotes::package_deps("riskParityPortfolio")

