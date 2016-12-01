\name{recountCells}
\alias{recountCells}

\title{Recount cells in each group}
\description{Count the number of cells in hyperspheres across a specified marker space, nested within pre-defined groups of cells.}

\usage{
recountCells(x, markers, tol=0.5, filter=10L)
}

\arguments{
\item{x}{A cyData object containing cell assignments into hyperspheres, such as that produced by \code{\link{countCells}}.}
\item{markers}{A vector specifying the markers for which hyperspheres should be constructed.}
\item{tol}{A numeric scalar proportional to the hypersphere radius.}
\item{filter}{An integer scalar specifying the minimum count sum required to report a hypersphere.} 
}

\details{
Each row of \code{x} corresponds to an existing hypersphere across some high-dimensional space, to which a set of cells are assigned.
This function takes this set of cells and assigns them into hyperspheres in dimensions specified by \code{markers}.
Thus, these new hyperspheres are \dQuote{nested} within the existing hypersphere, because only the cells in the latter are assigned into the former.
The counts and cell assignments for the nested hyperspheres are then reported in an output cyData object.

This function is useful in situations where the markers have different purposes.
For example, \code{x} could be constructed using only markers that define cell type.
This yields groups of cells that correspond to putative cell types.
\code{recountCells} can then be applied with \code{markers} that define, e.g., the activation status within each cell type.
In general, it is most interesting to use \code{markers} that were \emph{not} used to construct \code{x}.
Otherwise, by definition, all the cells in each hypersphere would have similar marker intensities.

To speed up the search, this function will only construct nested hyperspheres that are centered on the same cells as the hyperspheres in \code{x}.
This reduces the amount of distance re-calculations that are required.
However, hyperspheres that are not in \code{x} will not be considered as centres of nested hyperspheres.
This is generally fine if \code{x} has only been subjected to simple filtering based on average abundance.
In such cases, this missing hyperspheres will be of low abundance and would be removed by \code{filter} anyway.

The \code{tol} and \code{filter} arguments are interpreted as described for \code{\link{countCells}}.
Note that the distance calculation uses only the number of markers in \code{markers}, not the full marker space used to construct \code{x}.
}

\value{
A cyData object containing counts and cell assignments for nested hyperspheres.
This follows the same format as the output from \code{\link{countCells}}, i.e., each row is a hypersphere and each column is a sample.
Some additional fields are also included:
\itemize{
\item An extra \code{group} field is present in the \code{rowData}, specifying the row of \code{x} corresponding to the group in which each hypersphere is nested.
\item A \code{reused} field is present in \code{markerData}, specifying which markers were used in the re-count.
\item An integer scalar \code{retol} is present in the \code{metadata}, containing the tolerances used in re-counting.
}
The median intensities per hypersphere in \code{intensities} are only reported for the markers in \code{markers} -- all other values are set to \code{NA}.
Any existing \code{metadata} or \code{rowData} are preserved, while \code{cellIntensities} and \code{cellData} are unchanged.
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

