% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/famscore_to_file_bdeu.R
\name{famscore_to_file_bdeu}
\alias{famscore_to_file_bdeu}
\title{Compute BDeu family score for categorical data.}
\usage{
famscore_to_file_bdeu(
  data,
  max_parent_size = ncol(data) - 1,
  file_out,
  ess = 1,
  nlev
)
}
\arguments{
\item{data}{an n-by-d data frame or data matrix.}

\item{max_parent_size}{limit on the maximum parent set size.}

\item{file_out}{name of the file to which the scores are written (".score"}

\item{ess}{equivalent sample size parameter.}
}
\value{
a file with the parent scores in GOBNILP format.
}
\description{
Computes BDeu family score for all possible parent-child combinations. The scores
are written to file in GOBNILP format.
}
