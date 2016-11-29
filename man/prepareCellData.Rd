\name{prepareCellData}
\alias{prepareCellData}

\title{Prepare mass cytometry data}

\description{Convert single-cell marker intensities from a mass cytometry experiment into a format for efficient counting.}

\usage{
prepareCellData(x, naive=FALSE, markers=NULL, ...)
}

\arguments{
\item{x}{A named list of numeric matrices, where each matrix corresponds to a sample and contains expression intensities for each cell (row) and each marker (column).
Alternatively, a ncdfFlowSet object containing the same information.}
\item{naive}{A logical scalar specifying whether k-means clustering should be performed.}
\item{markers}{A vector specifying the markers to use in distance calculations.}
\item{...}{Additional arguments to pass to \code{\link{kmeans}}.}
}

\details{
This function constructs a cyData object from the marker intensities of each cell in one or more samples.

If \code{naive=FALSE}, this function performs k-means clustering on all the cells based on their marker intensities.
The number of clusters is set to the square-root of the total number of cells.
The cluster centres and cell assignments are then stored for later use in speeding up high-dimensional searches.
Intensity matrices from several samples are also merged into a single matrix for greater efficiency.

Note that \code{naive} does \emph{not} change the results of downstream functions, only the computational algorithm with which they are obtained.

If \code{markers} is specified, only the specified markers will be used in the distance calculations.
This also applies to calculations in downstream functions like \code{\link{countCells}} and \code{\link{neighborDistances}}.
All other markers will be ignored unless their usage is explicitly requested.
By default, \code{markers=NULL} which means that all supplied markers will be used in the calculations.
}

\value{
A cyData object with marker intensities for each cell stored in the \code{cellIntensities} slot.
In addition:
\itemize{
\item Marker names are stored as the row names of the \code{markerData} slot.
An additional \code{used} field also specifies whether the marker was used or not.
\item Sample names are stored as the column names of the output object.
\item \code{cellData} contains \code{sample.id}, an integer vector specifying the element of \code{x} that each cell was taken from;
    and \code{cell.id}, an integer vector specifying the row index of each cell in its original matrix.
\item If \code{naive=FALSE}, the metadata contains \code{cluster.centers}, a numeric matrix containing the centre coordinates of each cluster (column) for each marker (row);
    and \code{cluster.info}, a list containing information for each cluster.
}
Each element of \code{cluster.info} is a list, containing the zero-indexed column index of the output matrix that specifies the first cell in the cluster;
    as well as a numeric vector of distances between each cell in the cluster and the cluster centre.
Cells in \code{cellIntensities} are arranged in blocks corresponding to the clusters and ordered such that the distances are increasing.
}

\author{
Aaron Lun
}

\examples{
### Mocking up some data: ###
nmarkers <- 40
marker.names <- paste0("X", seq_len(nmarkers))
nsamples <- 10
sample.names <- paste0("Y", seq_len(nsamples))

x <- list()
for (i in sample.names) {
    ex <- matrix(rgamma(nmarkers*1000, 2, 2), ncol=nmarkers, nrow=1000)
    colnames(ex) <- marker.names
    x[[i]] <- ex
}

### Running the function: ###
cd <- prepareCellData(x)
cd
}

\seealso{
\code{\link{countCells}},
\code{\link{neighborDistances}}
}
