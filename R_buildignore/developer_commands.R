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
## Downloads (https://ipub.com/dev-corner/apps/r-package-downloads/)
##
library(cranlogs)
library(ggplot2)
library(dplyr)
library(tidyr)
downloads <- cran_downloads(from = "2017-12-20", 
                            package = c("riskParityPortfolio", "sparseIndexTracking", "sparseEigen", "spectralGraphTopology"))
downloads <- downloads %>% group_by(package) %>% mutate("cum_count" = cumsum(count)) %>% ungroup()
downloads %>% 
  filter(package == "riskParityPortfolio" & date >= "2018-12-25") %>%
  select("date", "count", "cum_count") %>%
  tail(8)
#ggplot(downloads, aes(x = date, y = count, color = package)) + geom_line() + ggtitle("Downloads")
#ggplot(downloads, aes(x = date, y = cum_count, color = package)) + geom_line() + ggtitle("Cumulative downloads")
labels <- c(count = "Number of downloads", cum_count = "Cumulative number of downloads")
downloads %>%
  gather("count", "cum_count", key = "count_type", value = "value") %>%
  ggplot(aes(x = date, y = value, color = package)) + 
  geom_line() +
  facet_wrap(~ count_type, ncol = 1, scales = "free", labeller = labeller(count_type = labels)) +
  ggtitle("Downloads") + xlab(NULL) + ylab(NULL) #+ ggsave("downloads.pdf", device = "pdf", scale = 0.5)




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


# Code tests (https://codecov.io/gh/mirca/riskParityPortfolio)
devtools::test()
covr::package_coverage()  # coverage of tests
#goodpractice::gp()  # overall checks


# CRAN check and submission (http://r-pkgs.had.co.nz/release.html)
#  checklist: https://kalimu.github.io/post/checklist-for-r-package-submission-to-cran/
devtools::check()
rcmdcheck::rcmdcheck()
devtools::build()
#devtools::revdep(pkg = "riskParityPortfolio")  # to check reverse dependencies
#devtools::build_win()  #to check under windows
#R CMD build .  # this is to generate tarball
#R CMD build . --compact-vignettes=gs+qpdf  
#R CMD check riskParityPortfolio_0.1.2.9000.tar.gz --as-cran  # this is before submission to CRAN
#R CMD install riskParityPortfolio_0.1.2.9000.tar.gz
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html


## Reverse dependencies
# tools::dependsOnPkgs("riskParityPortfolio")
# tools::package_dependencies("riskParityPortfolio", reverse = TRUE)
# devtools::revdep(pkg = "riskParityPortfolio") 
# devtools::revdep_check()
# remotes::package_deps("riskParityPortfolio")

