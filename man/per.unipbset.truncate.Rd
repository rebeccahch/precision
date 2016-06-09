% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/per.unipbset.truncate.R
\name{per.unipbset.truncate}
\alias{per.unipbset.truncate}
\title{Classification analysis of uniformly-handled data}
\usage{
per.unipbset.truncate(data, pbset.id = NULL, num.per.unipbset = 10)
}
\arguments{
\item{data}{expression data, rows as probes, columns as samples.}

\item{pbset.id}{unique probe-set name; default is NULL, the rownames of the dataset.}

\item{num.per.unipbset}{number of probes for each unique probe-set; default is 10.}
}
\value{
benchmark analysis results with list of models built and internal and external misclassification error stored, also a list of assignment stored
}
\description{
Performs classification analysis on the uniformly-handled data by reassigning samples to training and test set in Qin et al. (see reference).
}
\examples{
r.data.pl.p5 <- per.unipbset.truncate(data = r.data.pl,
num.per.unipbset = 5)
}
\references{
http://clincancerres.aacrjournals.org/content/20/13/3371.long
}
\keyword{data.setup}
