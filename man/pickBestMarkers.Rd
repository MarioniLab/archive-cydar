\name{pickBestMarkers}
\alias{pickBestMarkers}

\title{Pick best markers}

\description{Pick the best markers that distinguish between cells in and outside of a subpopulation.}

\usage{
pickBestMarkers(cell.data, chosen, tol=0.5, downsample=10, p=0.05, naive=FALSE)
}

\arguments{
\item{cell.data}{A matrix of marker intensities per cell produced by \code{\link{prepareCellData}}.}
\item{chosen}{A character vector specifying the names of the hyperspheres corresponding to the subpopulation of interest.}
\item{tol}{A numeric scalar, equal to that used in \code{\link{countCells}} to obtain the hypersphere counts.} 
\item{downsample}{A numeric scalar specifying the cell downsampling interval.}
\item{p}{A numeric scalar defining the quantiles for gating.}
\item{naive}{A logical scalar specifying whether a naive search should be used.}
}

\details{
The subpopulation is defined by a user-supplied set of hyperspheres in \code{chosen}.
By default, the name of each hypersphere from \code{countCells} specifies a row index of \code{cell.data}, corresponding to the cell on which the hypersphere is centred.
Provided these row names have not been altered, users can use them directly as \code{chosen}.

Cells in \code{cell.data} are downsampled according to \code{downsample}.
The function then identifies all cells in the downsampled set that were counted into any of the hyperspheres specified by \code{chosen} at the tolerance \code{tol}.
Note that \code{tol} \emph{must} be the same as the value used in \code{\link{countCells}} to originally count the cells into hyperspheres.
We also strongly recommended that \code{downsample} also be set to the same value.
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

# Selecting one of the mocked-up subpopulations.
mock <- kmeans(cnt$coordinates, centers=2)
X <- prcomp(cnt$coordinates)
plot(X$x[,1], X$x[,2], col=ifelse(mock$cluster==1L, "red", "blue"))

# Using the row names for selection (which are the column indices of 'cd')
selected <- rownames(cnt$coordinates)[mock$cluster==1L]

# Using the same default 'tol' as countCells.
pickBestMarkers(cd, selected)
}
