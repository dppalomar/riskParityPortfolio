---
layout: default
title: Home
nav_order: 1
description: "*riskParityPortfolio* is a software tool focused on the design of risk parity portfolios using fast, accurate, state-of-the-art optimization methods."
permalink: /
---

# Fast and scalable design of risk parity portfolios
{: .fs-9 }

**riskParityPortfolio** is a software tool focused on the design of risk parity
portfolios using fast, accurate, state-of-the-art optimization methods.

{: .fs-6 .fw-300 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/dppalomar/riskParityPortfolio){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Getting started

### Dependencies

The R version of **riskParityPortfolio** is build on top of awesome R packages including **Rcpp**,
**RcppEigen**, **quadprog**, **alabama**, and **nloptr**. All these packages can be installed via CRAN.

The Python version depends on **numpy**.

### Installation

1. The **stable** R version can be installed via CRAN as
```bash
install.packages("riskParityPortfolio")
```

2. The **development** R version can be installed via GitHub as
```bash
devtools::install_github("dppalomar/riskParityPortfolio")
```
<small>You must have previously installed the **devtools** package.</small>

3. The **stable** Python version can be installed via **pip** as
```bash
$ pip install riskparityportfolio
```

4. The **development** Python version can be installed via GitHub as
```bash
$ git clone https://github.com/dppalomar/riskParityPortfolio
$ cd python
$ pip install -e .
```

### Tutorials

- [See the package vignette](http://mirca.github.io/riskParityPortfolio/_static/getting_started.html) for a
detailed description of the mathematical methods that are available in **riskParityPortfolio**.

---

## About the project

**riskParityPortfolio** is developed on [GitHub](http://github.com/dppalomar/riskParityPortfolio)
by [Ze Vinicius](http://mirca.github.io) and [Daniel Palomar](http://www.danielppalomar.com).

### License

**riskParityPortfolio** is distributed by an
[GPL 3.0 License](https://github.com/dppalomar/riskParityPortfolio/blob/master/LICENSE).

### Contributing

We welcome all sorts of contributions. Please feel free to open an issue to report a bug or discuss a feature request in [our GitHub repo](https://github.com/dppalomar/riskParityPortfolio).

### Citation

If this package has been useful to you in any way, give us a star on [GitHub](http://github.com/dppalomar/riskParityPortfolio) :)
Additionally, if you've used **riskParityPortfolio** on your research, please consider citing the following resources:

- J. V. de M. Cardoso and D. P. Palomar (2019). riskParityPortfolio:
  Design of Risk Parity Portfolios. R package version 0.1.1. [https://CRAN.R-project.org/package=riskParityPortfolio](https://CRAN.R-project.org/package=riskParityPortfolio)
- Y. Feng, and D. P. Palomar (2015). SCRIP: Successive Convex Optimization Methods for Risk Parity Portfolio Design.
  IEEE Trans. on Signal Processing, vol. 63, no. 19, pp. 5285-5300. [https://doi.org/10.1109/TSP.2015.2452219](https://doi.org/10.1109/TSP.2015.2452219)
- F. Spinu (2013). An Algorithm for Computing Risk Parity Weights.
  [https://dx.doi.org/10.2139/ssrn.2297383](https://dx.doi.org/10.2139/ssrn.2297383)
- T. Griveau-Billion, J. Richard, and T. Roncalli (2013). A fast algorithm for computing high-dimensional
  risk parity portfolios. [https://arxiv.org/pdf/1311.4057.pdf](https://arxiv.org/pdf/1311.4057.pdf)
