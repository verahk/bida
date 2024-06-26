% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpt_arrays_from_bn.R
\name{cpt_arrays_from_bn}
\alias{cpt_arrays_from_bn}
\title{List CPTs of bn network}
\usage{
cpt_arrays_from_bn(bn)
}
\arguments{
\item{bn}{a \code{bn.fit.dnet} object (see \link[bnlearn:bn.fit]{bnlearn::bn.fit})}
}
\value{
a list with cpts stored as arrays
}
\description{
The CPTs associated with each \code{node} in a \code{bn.fit.dnet} object is available
\code{node$prob}.
To compute exact queries with \link{cpquery_from_cpt_arrays},
the scope of the CPTs must be given by the \code{dimnames} attribute.
}
\details{
\itemize{
\item \code{cpt_arrays_from_bn}
}
}
\keyword{internal}
