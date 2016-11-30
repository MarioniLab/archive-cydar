\name{recountCells}
\alias{recountCells}

\title{Recount cells in each group}
\description{Count the number of cells in hyperspheres across a specified marker space, nested within pre-defined groups of cells.}

\usage{
recountCells(x, markers, tol=0.5, downsample=NULL, filter=10, naive=FALSE)
}

\arguments{
\item{x}{A cyData object containing cell assignments into groups, such as that produced by \code{\link{countCells}}.}
\item{markers}{A vector specifying the markers for which hyperspheres should be constructed.}
\item{tol}{A numeric scalar proportional to the hypersphere radius.}
\item{downsample}{An integer scalar or vector specifying the frequency with which cells are sampled to form hyperspheres.}
\item{filter}{An integer scalar specifying the minimum count sum required to report a hypersphere.} 
\item{naive}{A logical scalar specifying whether a naive counting approach should be used.}
}

\details{
Each row of \code{x} corresponds to a pre-defined group of cells, usually those contained within a hypersphere when produced by \code{\link{countCells}}.
This function takes this group of cells and assigns them into hyperspheres in dimensions specified by \code{markers}.
Thus, these new hyperspheres are \dQuote{nested} within the pre-defined group, because only the cells in the latter are assigned into the former.
The counts and cell assignments for the nested hyperspheres are then reported in an output cyData object.

This function is useful in situations where the markers have different purposes.
For example, \code{x} could be constructed using only markers that define cell type.
This yields groups of cells that correspond to putative cell types.
\code{recountCells} can then be applied with \code{markers} that define, e.g., the activation status \emph{within} each cell type.

If \code{downsample=NULL}, an internal function is used to adjust the downsampling based on the number of cells within each group.
Groups with few cells get less downsampling so that the marker space is explored more thoroughly.
Groups with many cells get more downsampling, approaching \code{downsample=10}.
If a scalar or vector of length \code{nrow(x)} is specified, the supplied value will be used instead for each group.

The \code{tol}, \code{filter} and \code{naive} arguments are interpreted as described for \code{\link{countCells}}.
Note that the distance calculation uses only the number of markers in \code{markers}, not the full marker space used to construct \code{x}.
}

\value{
A cyData object containing counts and cell assignments for nested hyperspheres.
This follows the same format as the output from \code{\link{countCells}}, i.e., each row is a hypersphere and each column is a sample.
However, only information for the specified markers is stored in \code{intensities}, \code{cellIntensities} and \code{markerData}.
Also, an extra \code{group} field is present in the \code{rowData}, specifying the row of \code{x} corresponding to the group in which each hypersphere is nested.
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

