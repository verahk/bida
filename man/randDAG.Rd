% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randDAG.R
\name{randDAG}
\alias{randDAG}
\title{Draw a random adjacency matrix}
\usage{
randDAG(n, d, weighted = FALSE)
}
\arguments{
\item{n}{(integer) number of nodes}

\item{d}{(integer) expected number of neighbours}

\item{weighted}{(boolean) wheter each edge should be assigned weights. Default is FALSE.}
}
\value{
a n-by-n adjacency matrix
}
\description{
Wrapper around \link[pcalg:randDAG]{pcalg::randDAG} returning a adjacency matrix rather than a graphNEL object
}
\keyword{internal}
