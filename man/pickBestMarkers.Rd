\name{pickBestMarkers}
\alias{pickBestMarkers}

\title{Pick best markers}

\description{Pick the best markers that distinguish between cells in and outside of a set of hyperspheres.}

\usage{
pickBestMarkers(x, chosen, downsample=10, p=0.05, naive=FALSE)
}

\arguments{
\item{x}{A cyData object, constructed using \code{\link{countCells}}.}
\item{chosen}{A vector specifying the rows of \code{x} corresponding to the hyperspheres of interest.} 
\item{downsample}{A numeric scalar specifying the cell downsampling interval.}
\item{p}{A numeric scalar defining the quantiles for gating.}
\item{naive}{A logical scalar specifying whether a naive search should be used.}
}

\details{
A putative subpopulation is defined by a user-supplied set of hyperspheres in \code{chosen}.
Cells in \code{cellIntensities(x)} are downsampled according to \code{downsample}.
Then, this function identifies all cells in the downsampled set that were counted into any of the hyperspheres specified by \code{chosen} at the tolerance \code{tol}.
We recommend that \code{downsample} also be set to the same value as that used in \code{\link{countCells}} to construct \code{x}.
(This ensures that the identified cells are consistent with those that were originally counted.
It also avoids situations where no cells are counted into hyperspheres for rare subpopulations, which prevents GLM fitting as the response will only have one level.)

Relevant markers are identified by fitting a binomial GLM with LASSO regression to the downsampled cells, using the \code{\link{glmnet}} function.
The response is whether or not the cell was counted into the hyperspheres (and thus, the subpopulation).
The covariates are the marker intensities of each cell, used in a simple additive model with an intercept.
Upon fitting, the markers can be ranked from most to least important in terms of their ability to separate counted from uncounted cells.
This is done based on the LASSO iteration at which each marker's coefficient becomes non-zero - smaller values indicate more importance, while equal values indicate tied importance.
A panel of useful markers can subsequently be constructed by taking the top set from this ranking.

To evaluate the performance of each extra marker, we consider a progressive gating scheme.
For each marker, we define the gating boundaries as the interval between the \code{p} and \code{1-p} quantiles.
For a top set of markers, we calculate the number of cells from the subpopulation that fall inside the gating boundaries for each marker (i.e., true positives).
We repeat this for the number of cells not in the subpopulation (false positives).
This allows us to compute the recovery (i.e., sensitivity) of the gating scheme as the proportion of true positives out of the total number of cells in the subpopulation;
    and the contamination (i.e., non-specificity), as the proportion of false positives out of the total number of gated cells.
}

\value{
A data frame is returned, where each row is a marker ordered in terms of decreasing importance.
The contamination and recovery proportions of the top \code{n} markers are reported at row \code{n}, along with the LASSO iteration to denote ties.
The lower and upper gating boundaries are also reported for each marker.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{countCells}},
\code{\link{prepareCellData}},
\code{\link{glmnet}}
}

\examples{
# Mocking up some data with two clear subpopulations.
nmarkers <- 10L
ex1 <- matrix(rgamma(nmarkers*1000, 2, 2), ncol=nmarkers, nrow=1000)
ex2 <- ex1; ex2[,1:4] <- ex2[,1:4] + 1
ex <- rbind(ex1, ex2)
colnames(ex) <- paste0("X", seq_len(nmarkers))
cd <- prepareCellData(list(A=ex))
cnt <- countCells(cd, filter=1L)

# Selecting the first mocked-up subpopulations.
selected <- rowData(cnt)$center.cell <= 1000
pickBestMarkers(cnt, selected)
}
