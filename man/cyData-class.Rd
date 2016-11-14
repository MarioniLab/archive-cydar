\name{cyData-class}
\docType{class}

\alias{cyData}
\alias{cyData-class}
\alias{show,cyData-method}

\title{cyData class and methods}
\description{An overview of the cyData class and applicable methods.}

\usage{
cyData(markerData, medianIntensities=NULL, cellIntensities=NULL, assays=NULL, ...)
}

\arguments{
\item{markerData}{A DataFrame where each row corresponds to a marker, and is named according to that marker.}
\item{medianIntensities}{A numeric matrix of median intensities for each group of cells (row) and each marker (column).}
\item{cellIntensities}{A numeric matrix of intensities for each marker (row) and in each cell (column).}
\item{assays, ...}{Arguments to be passed to the \code{SummarizedExperiment} constructor.}
}

\details{
The cyData class is designed to store the cell counts for each group of cells (e.g., hyperspheres, clusters), along with the median intensities for each group.
It inherits from the SummarizedExperiment class and contains the additional slots:
\describe{
\item{markerData}{A DataFrame containing information about each marker in each row.
The row names are used as the marker names.
}
\item{medianIntensities}{A numeric matrix where each row corresponds to a group of cells and each column corresponds to a marker.
Each matrix entry contains the median intensity across all cells in a group for a particular marker.
}
\item{cellIntensities}{A numeric matrix containing intensities for each marker in each individual cell in the experiment.
Each column represents a cell while each row represents a marker.
This column-major setup is more amenable to fast processing later.
}
}
The above constructor will set \code{medianIntensities} and \code{assays} to a matrix with no rows, if not specified.
It will also set \code{cellIntensities} to a matrix with no columns and \code{cellData} to an empty list by default.
}

\value{
A cyData object containing the specified information.
}

\examples{
# A minimal example.
my.markers <- DataFrame(row.names=LETTERS)
cyd.minimal <- cyData(markerData=my.markers)

# Adding extra detail.
counts <- matrix(rpois(1000, 10), ncol=10)
medians <- matrix(rgamma(1000, 1, 1), ncol=10)
cell.int <- matrix(rgamma(10000, 1, 1), nrow=10)
marker.data <- DataFrame(row.names=LETTERS[1:10])
cyd <- cyData(counts, markerData=marker.data,
    medianIntensities=medians, cellIntensities=cell.int)
}

\author{
Aaron Lun
}
