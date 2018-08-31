<!-- README.md is generated from README.Rmd. Please edit that file -->
sparseIndexTracking
===================

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/sparseEigen)](http://cran.r-project.org/package=sparseEigen) [![CRAN Downloads](http://cranlogs.r-pkg.org/badges/sparseEigen)](http://cran.r-project.org/package=sparseEigen) ![CRAN Downloads Total](http://cranlogs.r-pkg.org/badges/grand-total/sparseEigen?color=brightgreen)

Computation of sparse portfolios for financial index tracking, i.e., joint selection of a subset of the assets that compose the index and computation of their relative weights (capital allocation). The level of sparsity of the portfolios, i.e., the number of selected assets, is controlled through a regularization parameter. Different tracking measures are available, namely, the empirical tracking error (ETE), downside risk (DR), Huber empirical tracking error (HETE), and Huber downside risk (HDR). See vignette for a detailed documentation and comparison, with several illustrative examples.

The package is based on the paper:

K. Benidis, Y. Feng, and D. P. Palomar, "Sparse Portfolios for High-Dimensional Financial Index Tracking," *IEEE Trans. on Signal Processing*, vol. 66, no. 1, pp. 155-170, Jan. 2018. (<https://doi.org/10.1109/TSP.2017.2762286>)

Installation
------------

``` r
# Installation from CRAN
install.packages("sparseIndexTracking")

# Installation from GitHub
# install.packages("devtools")
devtools::install_github("dppalomar/sparseIndexTracking")

# Getting help
library(sparseIndexTracking)
help(package = "sparseIndexTracking")
package?sparseIndexTracking
?spIndexTrack

# Citing this work
citation("sparseIndexTracking")
```

Vignette
--------

For more detailed information, please check the vignette: [GitHub-html-vignette](https://rawgit.com/dppalomar/sparseIndexTracking/master/vignettes/SparseIndexTracking-vignette.html), [GitHub-pdf-vignette](https://rawgit.com/dppalomar/sparseIndexTracking/master/vignettes/SparseIndexTracking-vignette.pdf), [CRAN-pdf-vignette](https://cran.r-project.org/web/packages/sparseIndexTracking/vignettes/SparseIndexTracking-vignette.pdf).

Usage of `spIndexTrack()`
-------------------------

Links
-----

Package: [CRAN](https://cran.r-project.org/package=sparseIndexTracking) and [GitHub](https://github.com/dppalomar/sparseIndexTracking).

README file: [GitHub-readme](https://rawgit.com/dppalomar/sparseIndexTracking/master/README.html) and [CRAN-readme](https://cran.r-project.org/web/packages/sparseIndexTracking/README.html).

Vignette: [GitHub-html-vignette](https://rawgit.com/dppalomar/sparseIndexTracking/master/vignettes/SparseIndexTracking-vignette.html) and [GitHub-pdf-vignette](https://rawgit.com/dppalomar/sparseIndexTracking/master/vignettes/SparseIndexTracking-vignette.pdf), [CRAN-pdf-vignette](https://cran.r-project.org/web/packages/sparseIndexTracking/vignettes/SparseIndexTracking-vignette.pdf).
