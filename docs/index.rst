.. riskParityPortfolio documentation master file, created by
   sphinx-quickstart on Sat Nov 10 08:44:52 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation
=============

``riskParityPortfolio`` is an `R` package for designing portfolios
that meet the risk parity criteria. It provides the implementation of the following
algorithms:
1) a Newton method for simple design of risk parity portfolios proposed by
`Spinu (2013) <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2297383>`_,
2) a variant of the cyclical coordinate descent algorithm for high dimensional portfolios proposed by
`Griveau-Billion (2013) <https://arxiv.org/pdf/1311.4057.pdf>`_ and applied to the risk parity formulation
presented in `Spinu (2013) <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2297383>`_,
3) general purpose non-linear constrained optimization solvers such as ``alabama`` and ``slsqp``, and
4) the successive convex approximation (SCA) proposed by
`Feng & Palomar (2016) <http://www.ece.ust.hk/~palomar/Publications_files/2015/FengPalomar-TSP2015%20-%20risk_parity_portfolio.pdf>`_.

Please, see the `Getting started <_static/getting_started.html>`_
tutorial for an introduction to risk parity portfolio design in R.

Installation
------------

The *development* version can be installed from GitHub as follows:

.. highlight:: r

::

   > devtools::install_github("dppalomar/riskParityPortfolio")

Citation
--------

If you made use of this package on your research, please consider citing the following resources:

- J. V. de M. Cardoso and D. P. Palomar (2018). riskParityPortfolio: Design of
  Risk Parity Portfolios. R package version 0.1.0. https://CRAN.R-project.org/package=riskParityPortfolio.
- Y. Feng, and D. P. Palomar, "SCRIP: Successive Convex Optimization Methods for
  Risk Parity Portfolio Design". IEEE Trans. on Signal Processing, vol. 63, no. 19,
  pp. 5285-5300, Oct. 2015. doi:10.1109/TSP.2015.2452219.
- F. Spinu, "An Algorithm for Computing Risk Parity Weights" (July 30, 2013).
  Available at SSRN: http://dx.doi.org/10.2139/ssrn.2297383.
- Griveau-Billion T., Richard J., and Roncalli T. "A fast algorithm for computing
  High-dimensional risk parity portfolios", ArXiv preprint: https://arxiv.org/pdf/1311.4057.pdf,
  2013.

Bug Reports
-----------

If you found a bug, please consider opening an issue ticket at our GitHub `repo <https://github.com/dppalomar/riskParityPortfolio/issues>`_.
