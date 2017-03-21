\name{CyData-subset}

\alias{[,CyData,ANY-method}
\alias{[,CyData,ANY,ANY-method}
\alias{[,CyData,ANY,ANY,ANY-method}
\alias{[<-,CyData,ANY,ANY,CyData-method}
\alias{subset,CyData-method}

\alias{cbind,CyData-method}
\alias{rbind,CyData-method}
\alias{c,CyData-method}

\title{Subsetting and combining CyData}
\description{Methods to subset and combine objects of the CyData class.}

\usage{
\S4method{[}{CyData,ANY,ANY}(x, i, j, ..., drop=TRUE)
\S4method{[}{CyData,ANY,ANY,CyData}(x, i, j) <- value
\S4method{subset}{CyData}(x, i, j)

\S4method{rbind}{CyData}(..., deparse.level=1)
\S4method{c}{CyData}(x, ..., recursive=FALSE)
\S4method{cbind}{CyData}(..., deparse.level=1)
}

\arguments{
\item{x}{A CyData object.}
\item{i, j}{
    A vector of subscripts, indicating the rows and columns to be subsetted for \code{i} and \code{j}, respectively.
    Rows correspond to cell-groups while columns correspond to samples.
}
\item{...}{
    For \code{cbind}, \code{...} contains CyData objects to be combined column-wise.
    For \code{rbind} and \code{c}, \code{...} contains CyData objects to be combined row-wise.
    For \code{c}, any objects are additional to that already specified in \code{x}.
    For other methods, this argument is ignored.
}
\item{drop}{
    A logical scalar, ignored.
}
\item{value}{
    A CyData object to replace the entries in the specified \code{i} and \code{j}.
}
\item{deparse.level, recursive}{Addition arguments, ignored.}
}

\details{
Subsetting of CyData objects consider groups of cells (e.g., clusters, hyperspheres) as the rows, and samples as the columns.
Only the rows of \code{intensities} are affected during row-wise subsetting.
Marker and cell information in \code{markerData}, \code{cellData} and \code{cellIntensities} is not modified.

A similar principle applies when combining different CyData objects, i.e., rows are groups of cells and columns are samples.
Values of \code{markerData} and \code{cellIntensities} should be identical in the objects to be combined.
Furthermore, for \code{cbind}, values of \code{intensities} should also be identical.
Note that \code{c} is a synonym for \code{rbind}.

When replacing entries with \code{[<-}, \code{markerData}, \code{cellData} and \code{cellIntensities} should be identical in \code{x} and \code{value}.
Column replacement also requires \code{intensities} and \code{cellAssignments} to be the same.
Simultaneous row and column replacement is not allowed as this will introduce inconsistencies into the sample ID annotation.

% Specifically, a column might be referring to one sample for some rows, and another sample for other rows.
}

\value{
A CyData object with different rows or columns, depending on whether subsetting or merging was performed.
}

\author{
Aaron Lun
}

\examples{
example(CyData) # Mocking up an object.

cyd[1:5,]
cyd[,6:10]
rbind(cyd, cyd)
cbind(cyd, cyd)
}
