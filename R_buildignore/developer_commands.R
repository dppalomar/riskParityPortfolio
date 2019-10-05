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



##
## Developer commands (https://r-pkgs.org/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::install()
#devtools::install(build_vignettes = TRUE)
library(riskParityPortfolio)

# Documentation
devtools::document()  #or Ctrl-Shift-L (to generate all documentation via roxygen)
?riskParityPortfolio


# Code tests (https://codecov.io/gh/mirca/riskParityPortfolio)
devtools::test()
covr::package_coverage()  # coverage of tests


# Compress pdf vignettes
#tools::compactPDF("vignettes/slides-RFinance2019.pdf", gs_quality = "ebook")
#tools::compactPDF("vignettes/slides-ConvexOptimizationCourseHKUST.pdf", gs_quality = "ebook")


# CRAN check and submission (http://r-pkgs.had.co.nz/release.html)
#  checklist: https://kalimu.github.io/post/checklist-for-r-package-submission-to-cran/
devtools::check()
rcmdcheck::rcmdcheck()
devtools::build()
#devtools::revdep(pkg = "riskParityPortfolio")  # to check reverse dependencies
#devtools::check_win_release()  #to check under windows
#R CMD build .  # this is to generate tarball
#R CMD build . --compact-vignettes=gs+qpdf  
#R CMD check riskParityPortfolio_0.2.1.tar.gz --as-cran --run-donttest  # this is before submission to CRAN
#R CMD install riskParityPortfolio_0.2.1.tar.gz
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html


## Reverse dependencies
# tools::dependsOnPkgs("riskParityPortfolio")
# tools::package_dependencies("riskParityPortfolio", reverse = TRUE)
# devtools::revdep(pkg = "riskParityPortfolio") 
# revdepcheck::revdep_check()  #source("https://install-github.me/r-lib/revdepcheck")
# remotes::package_deps("riskParityPortfolio")

