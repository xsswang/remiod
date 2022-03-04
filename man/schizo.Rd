
\name{schizo}

\alias{schizo}

\docType{data}

\title{National Institute of Mental Health shizophrenia study}

\description{Schizophrenia data from a randomized controlled trial with
  patients assigned to either drug or placebo group. "Severity of
  Illness" was measured, at weeks 0, 1, \ldots, 6, on a four category
  ordered scale. Most of the observations where made on weeks 0, 1, 3,
  and 6.}

\usage{data(schizo)}

\format{
  A data frame with 1603 observations on 437 subjects. Five vectors
  contain information on
  \describe{
    \item{\code{id}}{patient ID.}
    \item{\code{imps79}}{original response measurements on a numerical
      scale.}
    \item{\code{imps79b}}{binary response based on the cut-off value
      of 3.5 to the measurements on a numerical scale: 0 = normal to
      mildly ill and 1 = moderately to extremely ill}
    \item{\code{imps79o}}{ordinal response on a 4 category scale,
      "normal or borderline mentally ill" < "mildly or moderately ill",
      "markedly ill", "severely or among the most extremely ill".}
    \item{\code{tx}}{treatment indicator: 1 for drug, 0 for placebo.}
    \item{\code{week}}{week.}
  }
}

\details{The documentation file was copied from the \pkg{mixcat} package
  and slightly modified.}

\source{\url{https://hedeker.people.uic.edu/ml.html}}

\references{Hedeker, D. and R. Gibbons (2006). \emph{Longitudinal Data
  Analysis}. New Jersey, USA: John Wiley & Sons.}

\keyword{datasets}
