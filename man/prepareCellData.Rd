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
A numeric matrix where rows represent markers and columns represent cells (but not necessarily in the same order as the input cells in \code{x}).
Contiguous column-wise chunks of this matrix contain cells from the same cluster.
This matrix also contains several attributes:
\describe{
\item{\code{markers}:}{A character vector containing the marker names.}
\item{\code{samples}:}{A character vector containing the sample names.}
\item{\code{sample.id}:}{An integer vector specifying the sample matrix of \code{x} that each cell was taken from (zero indexed).}
\item{\code{cell.id}:}{An integer vector specifying the index of each cell in its original matrix (zero indexed).}
\item{\code{cluster.centers}:}{A numeric vector containing the centre coordinates of each cluster (column) for each marker (row).}
\item{\code{cluster.info}:}{A list containing information for each cluster.
Each element is a list, containing the zero-indexed column index of the output matrix that specifies the first cell in the cluster;
    and a numeric vector of distances between each cell in the cluster and the cluster centre.
Cells in the output matrix are ordered to correspond to the distances, in a manner such that the distances are increasing.
}
}
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
dim(cd)
names(attributes(cd))
}

\seealso{
\code{\link{countCells}}
}
