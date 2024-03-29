% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ladle.R
\name{ladle}
\alias{ladle}
\title{ladle}
\usage{
ladle(
  x = x,
  y = y,
  nslices = nslices,
  nboot = nboot,
  method = method,
  order = order
)
}
\arguments{
\item{x, }{predictors}

\item{y, }{response}

\item{nslices, }{parameter needed in the SIR methods}

\item{nboot, }{bootstrap 次数}

\item{method, }{AFM, pca, cca, ica, sir, save, dr, DEE, kernel sir.}

\item{order, }{parameter needed in the SIR methods}
}
\description{
The bootstrap ladle method by Ye
}
