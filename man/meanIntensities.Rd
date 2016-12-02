\name{meanIntensities}
\alias{meanIntensities}

\title{Compute mean marker intensities}
\description{Calcalute the mean intensity across cells in each group and sample for the specified markers.}

\usage{
meanIntensities(x, markers) 
}

\arguments{
\item{x}{A cyData object where each row corresponds to a group of cells, such as that produced by \code{\link{countCells}}.}
\item{markers}{A vector specifying the markers for which mean intensities should be calculated.}
}

\value{
A cyData object is returned equivalent to \code{x}, but with numeric matrices of mean intensities as additional elements of the \code{Assays} slot.
}

\details{
For each group of cells, the mean intensity across all assigned cells in each sample is computed.
This is returned as a matrix of mean intensities, with one value per sample (column) and hypersphere (row).
If a sample has no cells in a group, the corresponding entry of the matrix will be set to \code{NA}.

The groups in \code{x} should be defined using a different set of markers than in \code{markers}.
If the same markers were used for both functions, then a shift is unlikely to be observed.
This is because, by definition, the groups will contain cells with similar intensities for the markers used.

The idea is to use these values for weighted linear regression to identify a shift in intensity within each hypersphere.
The assumption is that the means are normally distributed, which is reasonable following the central limit theorem.
Means computed from more cells are also more precise, so the count for each group/sample can be used as a precision weight.
See the Examples for how this data can be used for entry into analysis packages like limma.

The mean intensity is used rather than the median, as the former is more sensitive to changes in intensity for a subset of cells assigned to a group.
Protection against random outliers is provided during linear regression, as these will increase the estimated variance.
This subsequently reduces significance of any shifts between conditions.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{countCells}}
}

\examples{
### Mocking up some data: ###
nmarkers <- 21
marker.names <- paste0("X", seq_len(nmarkers))
nsamples <- 5
sample.names <- paste0("Y", seq_len(nsamples))

x <- list()
for (i in sample.names) {
    ex <- matrix(rgamma(nmarkers*1000, 2, 2), ncol=nmarkers, nrow=1000)
    colnames(ex) <- marker.names
    x[[i]] <- ex
}

### Processing it beforehand with one set of markers: ###
cd <- prepareCellData(x, markers=marker.names[1:10])
cnt <- countCells(cd, filter=5)

## Computing the mean intensity for one marker: ###
cnt2 <- cydar:::meanIntensities(cnt, markers=marker.names[21])
library(limma)
mean.int.21 <- assay(cnt2, "mean.X21") 
cell.count <- assay(cnt2, "counts") 
el <- new("EList", list(E=mean.int.21, weight=cell.count))
}
