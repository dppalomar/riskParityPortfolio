.. riskParityPortfolio documentation master file, created by
   sphinx-quickstart on Sat Nov 10 08:44:52 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*************
Documentation
*************

The package `riskParityPortfolio` provides tools to design risk-parity portfolios.
In its simplest form, we consider the convex formulation with a unique solution
proposed by `Spinu (2013) <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2297383>`_
and use a cyclical method inspired by `Griveau-Billion (2013) <https://arxiv.org/pdf/1311.4057.pdf>`_.
For more general formulations, which are usually nonconvex, we implement the successive convex
approximation method proposed by `Feng & Palomar (2016) <http://www.ece.ust.hk/~palomar/Publications_files/2015/FengPalomar-TSP2015%20-%20risk_parity_portfolio.pdf>`_.

Please, see the `Getting started <_static/getting_started.html>`_
tutorial for an introduction to risk parity portfolio design in R.

Note that the Python version only supports the design of vanilla risk parity portfolios.

Installation
============

R
-

The *stable* version can be installed from CRAN as follows:

.. highlight:: r

::

   > install.packages("riskParityPortfolio")


The *development* version can be installed from GitHub as follows:

.. highlight:: r

::

   > devtools::install_github("dppalomar/riskParityPortfolio")

Python
------

The *stable* version can be installed from `pip` as follows:

.. highlight:: bash

::

   pip install riskparityportfolio

The *development* version can be installed from GitHub as follows:

.. highlight:: bash

::

   git clone https://github.com/dppalomar/riskParityPortfolio
   cd python
   pip install -e .

Citation
--------

If you made use of this package on your research, please consider citing the following resources:

- J. V. de M. Cardoso and D. P. Palomar (2019). riskParityPortfolio:
  Design of Risk Parity Portfolios. R package version 0.1.1. https://CRAN.R-project.org/package=riskParityPortfolio
- Y. Feng, and D. P. Palomar (2015). SCRIP: Successive Convex Optimization Methods for Risk Parity Portfolio Design.
  IEEE Trans. on Signal Processing, vol. 63, no. 19, pp. 5285-5300. https://doi.org/10.1109/TSP.2015.2452219
- F. Spinu (2013). An Algorithm for Computing Risk Parity Weights. https://dx.doi.org/10.2139/ssrn.2297383
- T. Griveau-Billion, J. Richard, and T. Roncalli (2013). A fast algorithm for computing
  high-dimensional risk parity portfolios. https://arxiv.org/pdf/1311.4057.pdf

Bug Reports
-----------

If you found a bug, please consider opening an issue ticket at our GitHub
`repo <https://github.com/dppalomar/riskParityPortfolio/issues>`_.
