\name{intensityRanges}
\alias{intensityRanges}

\title{Define intensity ranges}
\description{Set the ranges of the marker intensities, to direct construction of the colour bar for plotting.}

\usage{
intensityRanges(x, p=0.01)
}

\arguments{
\item{x}{A CyData object produced by \code{\link{prepareCellData}}.}
\item{p}{A numeric scalar specifying the quantile at which intensities should be bounded.}
}

\details{
For each marker, intensities across all cells are used to calculate the \code{p} and \code{1-p} quantiles.
This defines the lower and upper bound, respectively, to use as the \code{irange} argument in \code{\link{plotCellIntensity}}.
The aim is to prevent extreme outliers from skewing the distribution of colours.
This would result in loss of resolution at non-outlier values.

Note that, while the bounds are defined at the quantiles \code{p} and \code{1-p}, the colour gradient will not be computed across the percentiles.
That is, the "middle" of the gradient will not represent the median cell intensity.
Rather, the colour gradient is computed from the lower and upper bounds, so the middle with represent the average of the bounds.
In short, users should label the colour bar with the bounded intensities, rather than with the values of \code{p} or \code{1-p}.
}

\value{
A matrix specifying the lower and upper bounds (rows) on the intensity for each marker (columns).
}

\author{
Aaron Lun
}

\seealso{
\code{\link{prepareCellData}},
\code{\link{plotCellIntensity}}
}

\examples{
example(prepareCellData, echo=FALSE) # Using the mocked-up data set.
bounds <- intensityRanges(cd)

# Plotting example (using a subset for speed).
cd.subset <- t(cellIntensities(cd)[,1:1000])
coords <- prcomp(cd.subset)
chosen.marker <- 5
plotCellIntensity(coords$x[,1], coords$x[,2], 
    intensity=cd.subset[chosen.marker,], 
    irange=bounds[,chosen.marker])
}
