\name{spatialFDR}
\alias{spatialFDR}

\title{Compute the spatial FDR}
\description{
Computed adjusted p-values for all hyperspheres, using a density-weighted version of the Benjamini-Hochberg method.
}

\usage{
spatialFDR(coords, pvalues, neighbors=50, bandwidth=NULL, naive=FALSE)
}

\arguments{
\item{coords}{A numeric matrix of hypersphere coordinates, containing the median intensity of each marker (column) in each hypersphere (row).}
\item{pvalues}{A numeric vector of p-values for each hypersphere.}
\item{neighbors}{An integer scalar specifying the number of neighbors with which to compute the bandwidth.}
\item{bandwidth}{A numeric scalar specifying the bandwidth for density estimation.}
\item{naive}{A logical scalar specifying whether a naive search should be performed.}
}

\details{
Consider the set of significant hyperspheres, distributed in some manner across the M-dimensional space (for M markers).
The aim is to control the FDR across the subspaces containing significant hyperspheres.
This is subtly different from controlling the FDR across the hypersphere themselves, which will skew the results for densely occupied subspaces. 

Control of the spatial FDR is achieved by weighting the hyperspheres inversely proportional to their local densities.
This downweights hyperspheres in dense subspaces while upweighting hyperspheres in sparse subspaces.
The computed weights are then used as frequency weights in the Benjamini-Hochberg method, to control the FDR across subspaces.

The local density is calculated using a tricube kernel and the specified \code{bandwidth}.
If unspecified, \code{bandwidth} is set to the median of the distances to the \code{neighbors}-closest neighbor for all hyperspheres.
This usually provides stable density estimates while maintaining sensitivity to fine-scale structure.

Setting \code{naive=TRUE} will perform a naive search for nearest neighbors, rather than the more efficient convex method.
However, this should only affect computational efficiency and should not change the final results.
}

\author{
Aaron Lun
}

\value{
A numeric vector of adjusted p-values for all hyperspheres.
}

\examples{
coords <- matrix(rgamma(10000, 2, 2), nrow=1000)
pvalues <- rbeta(nrow(coords), 1, 2)
out <- spatialFDR(coords, pvalues)
}
