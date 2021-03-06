% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.handling.effect.R
\name{estimate.handling.effect}
\alias{estimate.handling.effect}
\title{Estimated handling effects}
\usage{
estimate.handling.effect(uhdata, nuhdata)
}
\arguments{
\item{uhdata}{the uniformly-handled dataset.
The dataset must have rows as probes and columns as samples.}

\item{nuhdata}{the nonuniformly-handled dataset.
The dataset must have rows as probes and columns as samples and the same dimensions and
the same probe names as the uniformly-handled dataset.}
}
\value{
an estimation of the handling effects
}
\description{
Estimate handling effects of an array of the nonuniformly-handled dataset
by taking the difference between its data and the data of
its matched array in the uniformly-handled dataset.
}
\examples{
handling.effect <- estimate.handling.effect(uhdata = uhdata.pl, nuhdata = nuhdata.pl)

}
\keyword{data.setup}
