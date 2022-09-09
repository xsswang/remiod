
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <span style="color: blue;">remiod</span>: Reference-based Multiple Imputation of Longitudinal Binary and Ordinal Outcomes with non-ignorable missingness

<!-- badges: start -->

[![CRAN
Status](https://www.r-pkg.org/badges/version/remiod)](https://CRAN.R-project.org/package=remiod)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/remiod)](https://cran.r-project.org/package=remiod)
[![GPL-3.0](https://img.shields.io/github/license/xsswang/remiod?logo=GNU&logoColor=FFFFFF&style=flat-square)](https://cran.r-project.org/package=remiod)
[![R build
status](https://github.com/xsswang/remiod/workflows/R-CMD-check/badge.svg)](https://github.com/xsswang/remiod/actions)
<!-- badges: end -->

The package **remiod** provides functionality to perform controlled
multiple imputation of binary and ordinal response in the Bayesian
framework. Implemented are (generalized) linear regression models for
binary data and cumulative logistic models / ordered probit models for ordered categorical data
(Wang and Liu 2022). It is also possible to fit multiple models of mixed
types simultaneously. Missing values in variables(if present) will be imputed
automatically.

**remiod** has two algorithmic backend. One is
[JAGS](https://mcmc-jags.sourceforge.io/), with which the function
performs some preprocessing of the data and creates a JAGS model, which
will then automatically be passed to
[JAGS](https://mcmc-jags.sourceforge.io/) with the help of the R package
[**rjags**](https://CRAN.R-project.org/package=rjags). The another is
based on the method proposed by Tang (Tang 2018).

Besides the main modelling functions, **remiod** also provides functions
to summarize and visualize results.

## Installation

**remiod** Can be from
[CRAN](https://cran.r-project.org/package=remiod):

``` r
install.packages("remiod")
```

Or, it can be installed from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("xsswang/remiod")
```

## Main functions

**remiod** provides the following main functions:

``` r
remiod                      #processing data and implementing MCMC sampling
extract_MIdata              #extract imputed data sets
miAnalyze                   #Perform analyses using imputed data and pool results
```

Currently, methods **remiod** implements include missing at random
(`MAR`), jump-to-reference (`J2R`), copy reference (`CR`), and delta
adjustment (`delta`). For `method = "delta"`, argument `delta` should
follow to specify a numerical values used in delta adjustment. These
methods can be requested through `extract_MIdata()`, and imputed
datasets can be analyzed using `miAnalyze()`.

Functions `summary()`, `coef()`, and `mcmcplot()` provide a summary of
the posterior distribution under MAR and its visualization.

## Minimal Example

``` r
data(schizow)

test = remiod(formula = y6 ~ tx + y0 + y1 + y3, data = schizow,
              trtvar = 'tx', algorithm = 'jags', method="MAR",
              ord_cov_dummy = FALSE, n.adapt = 10, n.chains = 1,
              n.iter = 100, thin = 2, warn = FALSE, seed = 1234)

extdt = extract_MIdata(object=test, method="J2R",mi.setting=NULL, M=10, minspace=2)
result = miAnalyze(y6 ~ y0 + tx, data = extdt, pool = TRUE)
```

## Support

For any help with regards to using the package or if you find a bug
please create a [GitHub
issue](https://github.com/xsswang/remiod/issues).

## Reference

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Erler2021" class="csl-entry">

Erler, NS, D Rizopoulos, and EMEH Lesaffre. 2021. “JointAI: Joint
Analysis and Imputation of Incomplete Data in R.” *Journal of
Statistical Software* 100 (20): 1–56.
<https://doi.org/10.18637/jss.v100.i20>.

</div>

<div id="ref-tang2017" class="csl-entry">

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
