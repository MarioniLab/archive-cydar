\name{recountCells}
\alias{recountCells}

\title{Recount cells in each group}
\description{Count the number of cells in hyperspheres across a specified marker space, nested within pre-defined groups of cells.}

\usage{
recountCells(x, markers, tol=0.5)
}

\arguments{
\item{x}{A cyData object containing cell assignments into hyperspheres, such as that produced by \code{\link{countCells}}.}
\item{markers}{A vector specifying the markers for which hyperspheres should be constructed.}
\item{tol}{A numeric scalar proportional to the hypersphere radius.}
}

\details{
Each row of \code{x} corresponds to an existing hypersphere across some high-dimensional space, to which a set of cells are assigned.
This function extends the hypersphere into the dimensions specified by \code{markers}.
Thus, each new hypersphere is \dQuote{nested} within the existing hypersphere in \code{x}.
Only the cells in the latter are assigned into the former
(though obviously, some cells will not be assigned if they are too distant from the centre in the new dimensions).

This function allows for fast recounting in situations where the markers have different purposes.
For example, \code{x} could be constructed using cell type-specific markers to define cell types.
\code{recountCells} can then be applied with \code{markers} that define, e.g., the activation status within each cell type.
In general, it is most interesting to use \code{markers} that were \emph{not} used to construct \code{x}.
Otherwise, by definition, all the cells in each hypersphere would have similar marker intensities.
}

\section{Note on the radius calculation}{
The output of this function is designed to be equivalent to directly running \code{\link{countCells}} with both new and old markers.
However, to speed up the counting, only the cells already assigned to each hypersphere in \code{x} are considered for re-counting.
This has some consequences for the results, as the radius scales with respect to \code{tol} and the number of markers.

Specifically, when new markers are specified in \code{markers}, the radius must increase to accommodate the increase in dimensions.
However, the cells were originally counted with a radius proportional to the (square root of the) old number of markers and \code{metadata(x)$tol}.
If the radius now increases, but only pre-assigned cells are used for re-counting, then there will be cells that are missed in the re-counts.

Thus, to preserve equivalence with \code{\link{countCells}} output, \code{tol} is decreased so that the radius does not change with new \code{markers}.
This shows up as a warning specifying the effective tolerance that was used in during re-counting.
Users can also avoid this problem by using a higher radius when constructing \code{x}, such that the radius calculated from the \code{tol} here will be smaller.
}

\value{
A cyData object containing counts and cell assignments for nested hyperspheres.
This follows the same format as the output from \code{\link{countCells}}, i.e., each row is a hypersphere and each column is a sample.
Some fields are modified:
\itemize{
\item \code{counts}, \code{intensities} and \code{cellAssignments} contain the relevant values for the nested hyperspheres.
\item The \code{used} field in \code{markerData} is set to all markers used in both the original and re-counting.
\item The \code{tol} value in the \code{metadata} is set to the (effective) tolerance used in re-counting.
}
}

\seealso{
\code{\link{countCells}}
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

### Processing it beforehand with one set of markers: ###
cd <- prepareCellData(x, markers=marker.names[1:20])
cnt <- countCells(cd, filter=5)

### Processing it afterwards with another set of markers: ###
rcnt <- recountCells(cnt, markers=marker.names[21:22])
rcnt
}

