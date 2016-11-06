\name{findFirstSphere}
\alias{findFirstSphere}

\title{Identifies the first non-redundant hyperspheres}
\description{Tests whether each hypersphere is not redundant to (i.e., lies more than a threshold distance away from) another hypersphere with a lower p-value.}

\usage{
findFirstSphere(coords, pvalues, threshold=1, naive=FALSE)
}

\arguments{
\item{coords}{A numeric matrix of hypersphere coordinates (median locations for all markers), where rows correspond to hyperspheres and columns correspond to markers.}
\item{pvalues}{A numeric vector of p-values, one for each hypersphere/row of \code{coords}.}
\item{threshold}{A numeric scalar specifying the maximum distance between the locations of two redundant hyperspheres.}
\item{naive}{A logical scalar specifying whether naive counting should be performed.}
}

\details{
This function iterates across the set of hyperspheres, typically ordered by decreasing significance.
It will tag a hypersphere as being redundant if its location lies within \code{threshold} of the location of a higher-ranking hypersphere in all dimensions.
In this manner, the set of all DA hyperspheres can be filtered down to a non-redundant subset that is easier to interpret.

Note that the criterion for redundancy mentioned above is equivalent to a Chebyshev distance, rather than Euclidean.
This is easier to interpret, especially given that the median intensity is defined separately for each marker.
Unlike in \code{\link{countCells}}, the threshold is not scaled by the number of markers because each hypersphere location is computed as an average across cells.
This means that there is generally no need to account for extra distance due to noise between cells.

The default \code{threshold} of unity assumes that the intensities have been transformed to or near a log10 scale.
It means that one hypersphere must vary from another by at least one log10-unit (i.e., a 10-fold change in intensity) in at least one marker to be considered non-redundant.
This avoids reporting many hyperspheres that differ from each other by relatively small, uninteresting shifts in intensity.
Greater resolution can be obtained by decreasing this value, e.g., to 0.5.

Note that setting \code{naive=TRUE} will only change the speed of the algorithm, not the results.

% Does the spatial FDR control the FDR across the reported non-redundant locations?
% In general, using spatialFDR() is fine; the non-redundant locations can be considered as centres of non-overlapping hypercubic partitions of side length "threshold".
% Thus, controlling the spatial FDR should be similar to controlling the FDR across these partitions.
%       Of course, there are some edge cases where this reasoning could fail, e.g., with diffuse DA subpopulations and tight non-DA subpopulations.
% If both subpopulations are detected, the spatial FDR may be low, because the former supplies a lot of true positive volume.
% However, the FDR across the non-redundant hyperspheres might be high if each subpopulation is represented by a single non-redundant hypersphere.
%       Ideally, we would like to control the FDR across the non-redundant hyperspheres directly.
% However, this is not possible because the definition of non-redundancy depends on the p-value.
% Applying the BH method to the non-redundant p-values would fail to account for the implicit multiple testing in the surrounding space.
}

\value{
A logical vector indicating whether each of the hyperspheres in \code{coords} is non-redundant.
}

\author{
    Aaron Lun
}

\examples{
# Mocking up some data.
coords <- matrix(rnorm(10000, 2, sd=0.3), nrow=1000)
pval <- runif(1000)

# Keep most significant non-redundant ("first") hyperspheres.
findFirstSphere(coords, pval)
}
