% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/support.R
\name{vl_ffv}
\alias{vl_ffv}
\title{Compute the force-velocity factor from Van Leuwen et al. (1990) using a modified Hill model}
\usage{
vl_ffv(V, k = 0.25, Vmax = 15)
}
\arguments{
\item{V}{numeric, he instantaneous contraction velocity (in muscle lengths, ML, per second).}

\item{k}{numeric, the constant determining the shape of the Hill curve.}

\item{Vmax}{numeric, the maximum contraction velocity (in ML per second).}
}
\value{
A single numeric value of the factor.
}
\description{
Compute the force-velocity factor from \insertCite{van1990function;textual}{JawCloseR}. Used internally by \code{close_jaw} to scale muscle-force output according to shortening velocity.
}
\details{
\loadmathjax This function returns the force-velocity factor given by a modified Hill  model \insertCite{hill1938heat}{JawCloseR} under both concentric \eqn{V>0} and eccentric \eqn{V<0} conditions published by Van Leuwen et al. 1990:

\mjseqn{ F_{FV}=\cfrac{V_{max}-V}{V_{max}+V/k}, \text{when } V>0}

\mjseqn{ F_{FV}=1.8-0.8 \times \cfrac{1+V/V_{max}}{1-7.56 \times V/(k \times V_{max})}, \text{when }  V<0}

The defualt value for \eqn{k} or 0.25 corresponds to vertebrate white muscle according to \insertCite{alexander2003principles;textual}{JawCloseR}.
}
\references{
\insertAllCited{}
}
