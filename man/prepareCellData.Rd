\name{prepareCellData}
\alias{prepareCellData}

\title{Prepare mass cytometry data}

\description{Convert single-cell marker intensities from a mass cytometry experiment into a format for efficient counting.}

\usage{
prepareCellData(x, ...)
}

\arguments{
\item{x}{A named list of numeric matrices, where each matrix corresponds to a sample and contains expression intensities for each cell (row) and each marker (column).
Alternatively, a ncdfFlowSet object containing the same information.}
\item{...}{Additional arguments to pass to \code{\link{kmeans}}.}
}

\details{
This function performs k-means clustering on all the cells based on their marker intensities.
The number of clusters is set to the square-root of the total number of cells.
The cluster centres and cell assignments are then stored for later use in speeding up high-dimensional searches.
Intensity matrices from several samples are also merged into a single matrix for greater efficiency.
}

\value{
A cyData object with marker intensities for each cell stored in the \code{cellIntensities} slot.
In addition:
\itemize{
\item Marker names are stored as the row names of the \code{markerData} slot.
\item Sample names are stored as the column names of the output object.
\item \code{cellData} contains \code{sample.id}, an integer vector specifying the sample matrix of \code{x} that each cell was taken from;
    and \code{cell.id}, an integer vector specifying the index of each cell in its original matrix.
\item The metadata contains \code{cluster.centers}, a numeric matrix containing the centre coordinates of each cluster (column) for each marker (row);
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
\code{\link{countCells}}
}
