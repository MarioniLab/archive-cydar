\name{medIntensities}
\alias{medIntensities}

\title{Compute median marker intensities}
\description{Calcalute the median intensity across cells in each group and sample for the specified markers.}

\usage{
medIntensities(x, markers) 
}

\arguments{
\item{x}{A cyData object where each row corresponds to a group of cells, such as that produced by \code{\link{countCells}}.}
\item{markers}{A vector specifying the markers for which median intensities should be calculated.}
}

\value{
A cyData object is returned equivalent to \code{x}, but with numeric matrices of sample-specific median intensities as additional elements of the \code{Assays} slot.
}

\details{
For each group of cells, the median intensity across all assigned cells in each sample is computed.
This is returned as a matrix of median intensities, with one value per sample (column) and hypersphere (row).
If a sample has no cells in a group, the corresponding entry of the matrix will be set to \code{NA}.

The groups in \code{x} should be defined using a different set of markers than in \code{markers}.
If the same markers were used for both functions, then a shift is unlikely to be observed.
This is because, by definition, the groups will contain cells with similar intensities for the markers used.

The idea is to use these values for weighted linear regression to identify a shift in intensity within each hypersphere.
The weight for each group/sample is defined as the number of cells, i.e., the \code{"counts"} assay in \code{x}.
This accounts for the precision with which the median is estimated, under certain assumptions.
See the Examples for how this data can be prepared for entry into analysis packages like limma.

The median intensity is used rather than the mean to ensure that shifts are interpreted correctly.
For example, mean shifts can be driven by strong changes in a subset of cells that are not representative of the majority of cells in the group.
This could lead to misinterpretation of the nature of the shift with respect to the group's overall identity.
}

%%% Assumptions of the precision weights %%%
% The use of the cell count as the precision weight makes use of the fact that the squared standard error of the median is inversely proportional to the number of observations.
% The exact relationship is only well-defined for normally distributed values, though the proportionality seems to hold for most distributions.
%
% # Should be robust to N:
% N <- 10
% a <- matrix(rexp(N*1000, 10), ncol=N)
% var(apply(a, 1, median)) * N
% b <- cbind(matrix(rnorm(N/2*1000, 0), ncol=N/2), matrix(rnorm(N/2*1000, 2), ncol=N/2))
% var(apply(b, 1, median)) * N
%
% An additional assumption of the weighting is that intensities for replicate samples are sampled from the same distribution.
% This might seem obvious, but it's actually less so because of biological variability that affects all cells in a sample.
% Finally, the only difference between conditions must involve a shift in the location of the distribution.
% Otherwise, if the distribution changes between conditions, the standard error will not follow the relationship above.
%
% All that being said, some weighting is probably better than nothing.

\section{Choosing between counting strategies}{
In situations where markers can be separated into two sets (e.g., cell type and signalling markers), there are two options for analysis.
The first is to define groups based on the \dQuote{primary} set of markers, then use \code{medIntensities} to identify shifts in each group for each of the \dQuote{secondary} markers.
This is the best approach for detecting increases or decreases in marker intensity that affect a majority of cells in each group.

The second approach is to recount cells into new groups using \code{\link{recountCells}} to focus on each secondary marker.
This provides more power to detect changes in marker intensity that only affect a subset of cells in each group.
Such changes are also easier to interpret as any correlation with respect to the primary markers for the affected subset can be studied.

The second approach is also more useful if one is interested in identifying cells with concomitant changes in multiple secondary markers.
Indeed, if we were interested in studying changes in all combinations of second markers, we would effectively revert to the obvious approach of just using all markers for counting.
However, this tends to be less effective for studying changes in a specific marker, due to the loss of precision with increased dimensionality.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{countCells}},
\code{\link{recountCells}}
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

## Computing the median intensity for one marker: ###
cnt2 <- medIntensities(cnt, markers=marker.names[21])
library(limma)
median.int.21 <- assay(cnt2, "med.X21") 
cell.count <- assay(cnt2, "counts") 
el <- new("EList", list(E=median.int.21, weight=cell.count))
}
