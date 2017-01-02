\name{poolCells}
\alias{poolCells}

\title{Pool cells for pre-processing}

\description{
Construct a flowFrame object by pooling cells from multiple (barcoded) samples, for use in common transformation and gating.
}

\usage{
poolCells(x, equalize=TRUE, n=NULL)
}

\arguments{
\item{x}{A named list of numeric matrices, where each matrix corresponds to a sample and contains expression intensities for each cell (row) and each marker (column).
Alternatively, a ncdfFlowSet object containing the same information.}
\item{equalize}{A logical scalar specifying whether the same number of cells should be taken from each sample for pooling.
If \code{FALSE}, all cells are used from all samples.}
\item{n}{A numeric scalar specifying the number of cells to be used from each sample if \code{equalize=TRUE}.
If \code{NULL}, this is set to the number of cells in the smallest sample.}
}

\details{
The idea is to use the pooled set of cells to estimate common parameters such as transformation values and gating thresholds.
Otherwise, if these parameters were estimated separately for each sample, they may distort the comparisons between samples.
This function is typically used to generate an object for use in \code{\link[flowCore]{estimateLogicle}} or in various gating functions like \code{\link{outlierGate}}.
This yields parameter values that can be applied to the full set of cells in the original \code{x} object.
}

\value{
A flowFrame object containing cells pooled from all samples.
}


\seealso{
\code{\link[flowCore]{flowFrame}},
\code{\link{outlierGate}},
\code{\link[flowCore]{estimateLogicle}}
}

\author{
Aaron Lun
}

\examples{
### Mocking up some data: ###
set.seed(100)
nmarkers <- 40
marker.names <- paste0("X", seq_len(nmarkers))
nsamples <- 10
sample.names <- paste0("Y", seq_len(nsamples))

x <- list()
for (i in sample.names) {
    ex <- matrix(rexp(nmarkers*1000, 0.01), ncol=nmarkers, nrow=1000)
    colnames(ex) <- marker.names
    x[[i]] <- ex
}

### Running the function: ###
ff <- poolCells(x)
ff

### Using for estimation: ###
library(flowCore)
trans <- estimateLogicle(ff, colnames(ff))
ff <- transform(ff, trans) # or, apply to original data.
}

