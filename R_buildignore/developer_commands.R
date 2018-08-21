##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type="source")
# Installation from GitHub
devtools::install_github("dppalomar/spectralGraphTopology")
# Installation from CRAN
install.packages("spectralGraphTopology")
# Getting help
library(spectralGraphTopology)
help(package="spectralGraphTopology")
package?spectralGraphTopology
?learnGraphTopology


##
## Developer commands (http://r-pkgs.had.co.nz/)
##
library(devtools)
devtools::load_all()  #or Ctrl-Shift-L
devtools::install()
#devtools::build()  # to generate the installation file

# Documentation
devtools::document()  #to generate all documentation via roxygen
?learnGraphTopology