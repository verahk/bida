% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interv_probs.R
\name{interv_probs}
\alias{interv_probs}
\alias{interv_probs_from_cpt_arrays}
\alias{interv_probs_from_bn}
\alias{interv_probs_x_exact}
\title{Compute intervention probabilities}
\usage{
interv_probs(obj, dmat, method = c("exact", "bn"))

interv_probs_from_cpt_arrays(cpts, method = c("exact", "bn"))

interv_probs_from_bn(bn, method = c("exact", "bn"))

interv_probs_x_exact(cpts, x, y, dmat)
}
\arguments{
\item{obj}{a object defining the distribution, see \code{cpts} or \code{bn}}

\item{method}{either \code{"exact"} or \code{"bn"}}

\item{cpts}{a list of cpt_arrays of length \code{n}, where \code{names(dimnames(cpts[[i]]))} gives
scope of the variable \code{i}.}

\item{bn}{a object of class \code{bnlearn::bn.fit.dnet} of length \code{n}}

\item{dag}{adjacency matrix of \code{obj}}
}
\value{
a \code{n}-by-\code{n} matrix with intervention probabilities

a list of \code{length(y)} with interventional cpts
}
\description{
Compute intervention probabilities

Compute exact intervention probabilities under interventions on \code{x}
}
\details{
These functions computes interventional CPTs for every cause-effect pair.
\itemize{
\item If \code{method == "exact"} the exact intervention probabilities is computed.
The \code{interv_probs_from_bn} function then first extract the cpt_arrays from
the bn-network.
\item If \code{method == "bn"}, the parameters are estimated using \link[bnlearn:cpquery]{bnlearn::cpdist}.
The \code{interv_probs_from_cpt_arrays} then first construct a \code{bn.fit.dnet} object.
\code{interv_probs} loops over every node, identifies the true descendants of
the node and compute the intervention distributions. For non-descendants,
the interventional CPTs are filled with the marginal distribution of the
pre-intervention distribution.
}
}
