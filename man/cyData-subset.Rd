\name{cyData-subset}

\alias{[,cyData,ANY-method}
\alias{[,cyData,ANY,ANY-method}
\alias{[,cyData,ANY,ANY,ANY-method}
\alias{[<-,cyData,ANY,ANY,cyData-method}
\alias{subset,cyData-method}

\alias{cbind,cyData-method}
\alias{rbind,cyData-method}
\alias{c,cyData-method}

\title{Subsetting and combining cyData}
\description{Methods to subset and combine objects of the cyData class.}

\usage{
\S4method{[}{cyData,ANY,ANY}(x, i, j, ..., drop=TRUE)
\S4method{[}{cyData,ANY,ANY,cyData}(x, i, j) <- value
\S4method{subset}{cyData}(x, i, j)

\S4method{rbind}{cyData}(..., deparse.level=1)
\S4method{c}{cyData}(x, ..., recursive=FALSE)
\S4method{cbind}{cyData}(..., deparse.level=1)
}

\arguments{
\item{x}{A cyData object.}
\item{i, j}{
    A vector of subscripts, indicating the rows and columns to be subsetted for \code{i} and \code{j}, respectively.
    Rows correspond to cell-groups while columns correspond to samples.
}
\item{...}{
    For \code{cbind}, \code{...} contains cyData objects to be combined column-wise.
    For \code{rbind} and \code{c}, \code{...} contains cyData objects to be combined row-wise.
    For \code{c}, any objects are additional to that already specified in \code{x}.
    For other methods, this argument is ignored.
}
\item{drop}{
    A logical scalar, ignored.
}
\item{value}{
    A cyData object to replace the entries in the specified \code{i} and \code{j}.
}
\item{deparse.level, recursive}{Addition arguments, ignored.}
}

\details{
Subsetting of cyData objects consider groups of cells (e.g., clusters, hyperspheres) as the rows, and samples as the columns.
Only the rows of \code{medianIntensities} are affected during row-wise subsetting.
Marker and cell information in \code{markerData}, \code{cellData} and \code{cellIntensities} is not modified.

A similar principle applies when combining different cyData objects, i.e., rows are groups of cells and columns are samples.
Values of \code{markerData} and \code{cellIntensities} should be identical in the objects to be combined.
Furthermore, for \code{cbind}, values of \code{medianIntensities} should also be identical.
Note that \code{c} is a synonym for \code{rbind}.

When replacing entries with \code{[<-}, \code{markerData}, \code{cellData} and \code{cellIntensities} should be identical in \code{x} and \code{value}.
}

\author{
Aaron Lun
}

\examples{
example(cyData) # Mocking up an object.

cyd[1:5,]
cyd[,6:10]
rbind(cyd, cyd)
cbind(cyd, cyd)
}
