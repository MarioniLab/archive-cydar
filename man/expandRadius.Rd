\name{expandRadius}
\alias{expandRadius}

\title{Expand the hypersphere radius}

\description{Expands the hypersphere radius to account for intensity shifting between non-barcoded samples.}

\usage{
expandRadius(x, design=NULL, markers=NULL, tol=0.5)
}

\arguments{
\item{x}{A CyData object produced by \code{\link{prepareCellData}}.}
\item{design}{A numeric matrix specifying the experimental design.}
\item{markers}{A vector specifying the markers to use.}
\item{tol}{A numeric scalar proportional to the hypersphere radius, see \code{\link{countCells}}.}
}

\details{
This function increases the hypersphere radius to account for random shifts in marker intensity between non-barcoded samples.
The required increase is estimated by taking the mean of all intensities for each marker in each sample;
computing the variance of the mean intensities across samples for each marker;
and taking the mean variance across all markers.
This is equivalent to the square of the extra distance between cells caused by intensity shifts between samples.

The estimated increase is added onto \code{tol}, and the returned value can be directly used in the \code{tol} argument of \code{\link{countCells}}.
This expands the hyperspheres to ensure that corresponding subpopulations in different samples are still counted together.
Otherwise, an intensity shift in one sample may move the cells in a subpopulation out of a hypersphere.
This will inflate the variability if it occurs between replicate samples, and introduce spurious differences if it occurs between samples in different conditions.

If \code{markers} is specified, only the indicated markers are used to calculate the mean variance.
Otherwise, the markers to use are extracted from \code{markerData(x)$used}, which is set by \code{\link{prepareCellData}}.
By default, all markers are used unless otherwise requested.
}

\value{
A numeric scalar specifying a modified \code{tol} to use in \code{\link{countCells}}.
}

\examples{
### Mocking up some data: ###
nmarkers <- 20
marker.names <- paste0("X", seq_len(nmarkers))
nsamples <- 8
sample.names <- paste0("Y", seq_len(nsamples))

x <- list()
for (i in sample.names) {
    ex <- matrix(rgamma(nmarkers*1000, 2, 2), ncol=nmarkers, nrow=1000)
    ex <- t(t(ex) + rnorm(nmarkers, 0, 0.25)) # Adding a shift per marker
    colnames(ex) <- marker.names
    x[[i]] <- ex
}

### Running the function: ###
cd <- prepareCellData(x)
expandRadius(cd)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{prepareCellData}},
\code{\link{countCells}}
}
