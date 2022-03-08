
<!-- README.md is generated from README.Rmd. Please edit that file -->

# remiod: Reference-based Controlled Multiple Imputation of Longitudinal Binary and Ordinal Outcomes with non-ignorable missingness

<!-- badges: start -->

[![GPL-3.0](https://img.shields.io/github/license/xsswang/remiod?logo=GNU&logoColor=FFFFFF&style=flat-square)](https://github.com/xsswang/remiod/main/LICENSE)
[![R build
status](https://github.com/xsswang/remiod/workflows/R-CMD-check/badge.svg)](https://github.com/xsswang/remiod/actions)
<!-- badges: end -->

The package **remiod** provides functionality to perform controlled
multiple imputation of binary and ordinal response in the Bayesian
framework. Implemented are (generalized) linear regression models for
binary data and cumulative logistic models for ordered categorical data
(Wang and Liu 2022). It is also possible to fit multiple models of mixed
types simultaneously. Missing values in (if present) will be imputed
automatically.

**remiod** has two algorithmic backends. One is
[JAGS](https://mcmc-jags.sourceforge.io/), with which the function
performs some preprocessing of the data and creates a JAGS model, which
will then automatically be passed to
[JAGS](https://mcmc-jags.sourceforge.io/) with the help of the R package
[**rjags**](https://CRAN.R-project.org/package=rjags). The another is
based on the method proposed by Tang (Tang 2018).

Besides the main modelling functions, **remiod** also provides functions
to summarize and visualize results.

## Installation

you can install **remiod** from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("xsswang/remiod")
```

## Main functions

**remiod** provides the following main functions:

``` r
remiod                      #processing data and implementing MCMC sampling
extract_MIdata              #extract imputed data sets
```

Currently, methods **remiod** implements include missing at random ,
jump-to-reference , copy reference , and delta adjustment . For ,
argument should follow to specify a numerical values used in delta
adjustment

Functions `summary()`, `coef()`, and `mcmcplot()` provide a summary of
the posterior distribution and its visualization.

## Minimal Example

``` r
data(schizow)

test = remiod(formula = y6 ~ tx + y0 + y1 + y3, data = schizow,
              trtvar = 'tx', algorithm = 'jags', method="MAR",
              ord_cov_dummy = FALSE, n.adapt = 10, n.chains = 1,
              n.iter = 100, thin = 2, warn = FALSE, seed = 1234)

extdt = extract_MIdata(object=test, method="J2R",mi.setting=NULL, M=10, minspace=2)
```

## Reference

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Erler2021" class="csl-entry">

Erler, NS, D Rizopoulos, and EMEH Lesaffre. 2021. “JointAI: Joint
Analysis and Imputation of Incomplete Data in R.” *Journal of
Statistical Software* 100 (20): 1–56.
<https://doi.org/10.18637/jss.v100.i20>.

</div>

<div id="ref-tang2018" class="csl-entry">

Tang, Y. 2018. “Controlled Pattern Imputation for Sensitivity Analysis
of Longitudinal Binary and Ordinal Outcomes with Nonignorable Dropout.”
*Statistics in Medicine* 37 (9): 1467–81.
<https://doi.org/10.1002/sim.7583>.

</div>

<div id="ref-wang2022" class="csl-entry">

Wang, T, and Y Liu. 2022. “Remiod: Reference-Based Controlled Multiple
Imputation of Longitudinal Binary and Ordinal Outcomes with
Non-Ignorable Missingness.” *arXiv* 2203.02771.
<https://arxiv.org/pdf/2203.02771>.

</div>

</div>
