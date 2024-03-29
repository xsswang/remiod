% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_subset.R
\name{get_subset}
\alias{get_subset}
\title{Extract specific parameters from MCMC samples}
\usage{
get_subset(object, subset, warn = TRUE, mess = TRUE)
}
\arguments{
\item{object}{an object of class \code{MCMC}.}

\item{subset}{subset of parameters (columns of the mcmc object) to be used. See
https://nerler.github.io/JointAI/articles/SelectingParameters.html
for key-words of subseting parameters. Besides, \code{selected_parms}
and \code{selected_vars} are new key-words for arbitrarily selecting
parameters.}

\item{warn}{logical, should warning messages be displayed?}

\item{mess}{logical, should messages be displayed?}
}
\description{
Extract specific parameters from MCMC samples
}
\examples{
\donttest{
data(schizow)

test = remiod(formula = y6 ~ tx + y0 + y1 + y3, data = schizow,
              trtvar = 'tx', algorithm = 'jags', method="MAR",
              ord_cov_dummy = FALSE, n.adapt = 10, n.chains = 1,
              n.iter = 0, thin = 1, warn = FALSE, seed = 1234)

pms = c("beta[2]","alpha[2]","alpha[6]","alpha[9]")
mcsub = get_subset(object = test$mc.mar, subset=c(selected_parms = list(pms)))

}
}
