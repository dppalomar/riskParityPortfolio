.. riskParityPortfolio documentation master file, created by
   sphinx-quickstart on Sat Nov 10 08:44:52 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to riskParityPortfolio!
===============================

``riskParityPortfolio`` allows users to design portfolios that meet the risk parity criteria.
We solve the underlying optimization problem using three approaches: 1) a Newton method
for simple risk parity design proposed by `Spinu (2013) <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2297383>`_,
2) general purpose non-linear constrained optimization solvers such as ``alabama`` and ``slsqp``,
and 3) the successive convex approximation (SCA) proposed by
`Feng & Palomar (2016) <http://www.ece.ust.hk/~palomar/Publications_files/2015/FengPalomar-TSP2015%20-%20risk_parity_portfolio.pdf>`_.

See the `Getting started <_static/getting_started.html>`_ tutorial for a comparison
between approaches.

Installation
------------

The *development* version can be installed from GitHub as follows:

.. highlight:: r

::

   > devtools::install_github("dppalomar/riskParityPortfolio")

Citation
--------

If you made use of this package on your research, please cite the following works:

- Y. Feng, and D. P. Palomar, "SCRIP: Successive Convex Optimization Methods for
  Risk Parity Portfolio Design". IEEE Trans. on Signal Processing, vol. 63, no. 19,
  pp. 5285-5300, Oct. 2015. doi:10.1109/TSP.2015.2452219
- F. Spinu, "An Algorithm for Computing Risk Parity Weights" (July 30, 2013).
  Available at SSRN: http://dx.doi.org/10.2139/ssrn.2297383

Bug Reports
-----------

If you found a bug, please consider opening an issue ticket at our GitHub `repo <https://github.com/dppalomar/riskParityPortfolio/issues>`_.
