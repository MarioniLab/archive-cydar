\name{findFirstSphere}
\alias{findFirstSphere}

\title{Identifies the first non-redundant hyperspheres}
\description{Tests whether each hypersphere is not redundant to (i.e., lies more than a threshold distance away from) another higher-ranking hypersphere.}

\usage{
findFirstSphere(coords, threshold=1, naive=FALSE)
}

\arguments{
\item{coords}{A numeric matrix of hypersphere coordinates (median locations for all markers), where rows correspond to hyperspheres and are sorted by some metric.}
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

% One question that arises is that, now that you have "non-redundant" hyperspheres, why can't you control the FDR across them directly?
% The answer is because the process of defining non-redundancy requires searching for the lowest p-values.
% You don't account for the implicit multiple testing in the surrounding space, which is likely to result in some liberalness.
% This is equivalent to partitioning the space, and then taking the lowest p-value in each partition as the representative for testing.
%
% On the flipside, the other question is if we controlled the spatial FDR, what is the FDR across the non-redundant hyperspheres?
% It should be very similar, as you now have one entry per hypercubic partition (of size equal to 2*threshold).
% The problem above doesn't apply because we're not using the non-redundant p-values for testing.
% To extend the analogy, it's like selecting the best p-value as a representative for examination, but not for testing.
%
% There might be some loss of control due to edge effects (as the partitions are always centered around the most significant local hypersphere).
% However, these should largely cancel out between the true and false positives, so the FDR shouldn't be too far off.
}

\value{
A logical vector indicating whether each of the hyperspheres in \code{coords} is non-redundant.
}

\author{
    Aaron Lun
}

\examples{
coords <- matrix(rnorm(10000, 2, sd=0.3), nrow=1000)
pval <- runif(1000)

# Sort by significance.
o <- order(pval)
coords <- coords[o,]

# Keep most significant non-redundant ("first") hyperspheres.
findFirstSphere(coords)
}
