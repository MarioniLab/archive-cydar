\name{Plot Cells}
\alias{plotCellLogFC}
\alias{plotCellIntensity}

\title{Plot cell or hypersphere data}
\description{Visualize cells or hyperspheres in low-dimensional space, coloured by marker intensities or log-fold changes.}

\usage{
plotCellLogFC(x, y, logFC, max.logFC=NULL, zero.col=0.8, 
    length.out=100, pch=16, ...)

plotCellIntensity(x, y, intensity, irange=NULL, 
    length.out=100, pch=16, ...)
}

\arguments{
\item{x, y}{A numeric vector of coordinates for each feature (i.e., cell or hypersphere).}
\item{logFC}{A numeric vector of log-fold changes for each feature.}
\item{max.logFC}{A numeric scalar specifying the maximum absolute log-fold change.}
\item{zero.col}{A numeric scalar between 0 and 1, specifying the greyscale intensity to represent a log-fold change of zero.}
\item{intensity}{A numeric vector specifying the marker intensities for each feature.}
\item{irange}{A numeric vector of length 2, specifying the upper and lower bound for the intensities.}
\item{length.out}{An integer scalar specifying the resolution of the colour bar.}
\item{pch, ...}{Additional arguments to pass to \code{plot}.}
}

\details{
\code{plotCellLogFC} will colour the points from blue (negative log-FC) to grey (zero log-FC) to red (positive log-FC).
The darkness of the grey colour is set with \code{zero.col}.
If \code{max.logFC} is not \code{NULL}, extreme values in \code{logFC} are winsorized to lie within \code{[-max.logFC, max.logFC]}.
This preserves the resolution for smaller log-fold changes.

\code{plotCellIntensity} will colour the points using the viridis colour scheme, i.e., purple (low intensity) to green (medium) to yellow (high).
If \code{irange} is not \code{NULL}, extreme values in \code{intensity} will be winsorized to lie within \code{irange}.
Like before, this preserves the resolution for smaller changes in intensity. 
Users should consider using \code{\link{intensityRanges}} to define appropriate values of \code{irange} for each marker.
}

\value{
A vector of colours of length \code{length.out} is returned, containing the colour gradient used for plotting.
The vector names contains the numeric values associated with each colour.
This can be used to construct a colour bar.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{viridis}},
\code{\link{intensityRanges}}
}

\examples{
# Making up some coordinates.
x <- rnorm(100)
y <- rnorm(100) 

# Log-FC plot and colour bar.
logFC <- rnorm(100)
out <- plotCellLogFC(x, y, logFC)
out <- plotCellLogFC(x, y, logFC, max.logFC=0.5)

plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
start.loc <- seq(-0.5, 0.5, length.out=length(out))
rect(-0.5, start.loc, 0.5, start.loc+diff(start.loc)[1], col=out, border=NA)
text(0, -0.5, pos=1, names(out)[1], cex=1.2)
text(0, 0.5,  pos=3, names(out)[length(out)], cex=1.2)
text(-0.6, 0, srt=90, "Log-FC", cex=1.2)

# Intensity plot and colour bar.
intensities <- rgamma(100, 2, 2)
out <- plotCellIntensity(x, y, intensities)
out <- plotCellIntensity(x, y, intensities, irange=c(0, 2))

plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
start.loc <- seq(-0.5, 0.5, length.out=length(out))
rect(-0.5, start.loc, 0.5, start.loc+diff(start.loc)[1], col=out, border=NA)
text(0, -0.5, pos=1, names(out)[1], cex=1.2)
text(0, 0.5,  pos=3, names(out)[length(out)], cex=1.2)
text(-0.6, 0, srt=90, "Intensity", cex=1.2)
}

