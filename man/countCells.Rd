\name{countCells}
\alias{countCells}

\title{Count cells in high-dimensional space}

\description{
Count the number of cells from each sample lying inside hyperspheres in high-dimensional space.
}

\usage{
countCells(x, tol=0.5, BPPARAM=SerialParam(), downsample=10, filter=10, naive=FALSE)
}

\arguments{
\item{x}{A CyData object produced by \code{\link{prepareCellData}}.}
\item{tol}{A numeric scalar proportional to the hypersphere radius.}
\item{BPPARAM}{A BiocParallelParam object specifying how parallelization is to be performed.}
\item{downsample}{An integer scalar specifying the frequency with which cells are sampled to form hyperspheres.}
\item{filter}{An integer scalar specifying the minimum count sum required to report a hypersphere.} 
\item{naive}{A logical scalar specifying whether a naive counting approach should be used.}
}

\details{
Consider that each cell defines a point in M-dimensional space (where M is the number of markers), based on its marker intensities.
This function constructs hyperspheres and counts the number of cells from each sample lying within each hypersphere.
In this manner, the distribution of cells across the space can be quantified.
For each hypersphere, cell counts for all samples are reported along with the median intensity across the counted cells for each marker.

Each hypersphere is centered on a cell to ensure that only occupied spaces are counted.
However, for high-density spaces, this can result in many redundant hyperspheres.
To reduce computational work, only a subset of cells are used to define hyperspheres.
The downsampling frequency is specified by \code{downsample}, e.g., only every 10th cell is used to make a hypersphere by default.

Each hypersphere also has a radius of \code{tol*sqrt(M)} (this relationship avoids loss of counts as M increases).
\code{tol} can be interpreted as the acceptable amount of deviation in the intensity of a single marker for a given subpopulation.
The default value of 0.5 means that, for any one marker, cells with +0.5 or -0.5 intensity will be counted into the same subpopulation.
This value is sensible as intensities are usually on a log-10 scale, such that a total of 10-fold variability in marker intensities is tolerated.

The coordinates are reported as (weighted) medians across all cells in each hypersphere.
This better reflects the location of the hypersphere if the cells are not distributed around the centre.
Each cell is weighted inversely proportional to the total number of cells in the corresponding sample.
This ensures that large samples do not dominate the median calculation.

All hyperspheres with count sums below \code{filter} are removed by default.
Such hyperspheres do not have enough counts (and thus, information) for downstream analyses.
Removing them reduces the amount of memory required to form the output matrix.

The default setting of \code{naive=FALSE} will use a method similar to that described by Samusik et al. (2016) to speed up the search.
Here, the distance between the hypersphere and cluster centers is first computed to decide which, if any, cells in the cluster should be considered for counting.
If \code{naive=TRUE}, a slower naive approach will be used where distances are computed between all pairs of cells.
This is generally only useful for testing.

If \code{markerData(x)$used} is not all \code{TRUE}, only the specified markers will be used in the distance calculations.
Median coordinates will not be calculated for any discarded markers, instead being reported as \code{NA} in the utput \code{intensities}.
Note that the \code{used} field should not be set directly -- \code{\link{prepareCellData}} must be run again to change the set of used markers.

Users can also increase speed by setting \code{BPPARAM} to use multiple cores.
Neither this or \code{naive} will not affect the final results, only the speed with which they are obtained.
}

\value{
A CyData object containing the following information:
\describe{
\item{\code{counts}}{An integer matrix of counts for each hypersphere (row) and sample (column) in the \code{Assays} slot.}
\item{\code{intensities}:}{A numeric matrix of median intensities for each hypersphere (row) and marker (column), stored in the \code{intensities} slot.}
\item{\code{cellAssignments}:}{A list of integer vectors specifying the cells belonging in each group.
This is stored in a compressed format, see \code{\link{packIndices}} for details.}
\item{\code{sample.id}:}{An integer vector containing an integer ID for each sample.
This corresponds to values in \code{cellData(x)$sample.id}, to ensure that they can be properly interpreted regardless of column subsetting or combining.}
\item{\code{totals}:}{An integer vector specifying the total number of cells in each sample, stored as a field in the \code{colData} slot.}
\item{\code{tol}:}{An integer scalar equal to \code{tol}, stored in the \code{metadata} slot.}
\item{\code{center.cell}:}{An integer vector specifying the column of \code{cellIntensities} containing the centre of each hypersphere, stored as a field in the \code{rowData} slot.}
}
Any existing \code{rowData}, \code{intensities} and \code{assays} are discarded.
All other slots are left unchanged from the values supplied in \code{x}.
}

\author{
Aaron Lun
}

\references{
Samusik N, Good Z, Spitzer MH et al. (2016).
Automated mapping of phenotype space with single-cell data.
\emph{Nat. Methods} 13:493-496
}

\seealso{
\code{\link{prepareCellData}},
\code{\link{packIndices}}
}
\examples{
example(prepareCellData, echo=FALSE)
downsample <- 10L
tol <- 0.5

cnt <- countCells(cd, filter=1, downsample=downsample, tol=tol)
cnt
}
