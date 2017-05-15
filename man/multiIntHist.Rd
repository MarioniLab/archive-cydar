\name{multiIntHist}
\alias{multiIntHist}

\title{multiIntHist}
\description{Generate intensity histograms from multiple batches.}

\usage{
multiIntHist(collected, cols=NULL, xlab="Intensity", ylab="Density",
    lwd=2, lty=1, pch=16, cex=2, ...) 
}

\arguments{
\item{collected}{A list of numeric vectors, where each vector contains intensities for a given marker from all cells of a single batch.}
\item{cols}{A vector of R colours of the same length as \code{collected}, to be used in colouring the histograms.}
\item{xlab, ylab}{Strings specifying the x- and y-axis labels.}
\item{lwd, lty}{Parameters for plotting the histogram traces.}
\item{pch, cex}{Parameters for plotting the frequency of zeroes.}
\item{...}{Other arguments to pass to \code{plot}.}
}

\details{
A histogram is constructed for the set of intensities from each batch, and the histogram outline is plotted with the specified parameters.
The frequency of intensities at zero (or negative values) is indicated with a single point at an intensity of zero.
This ensures that the number of events at zero and small non-zero intensities can be distinguished.

The process is repeated for all batches so that intensity distributions can be compared between batches.
If \code{cols=NULL}, the \code{rainbow} colour palette is automatically used to generate the colour for each batch. 
Some small jitter is added to the zero points so that they do not completely overlap each other.
}

\value{
Histogram traces representing the intensity distributions are produced on the current graphics device.
}

\author{Aaron Lun}

\seealso{
\code{\link{normalizeBatch}}
}

\examples{
multiIntHist(list(rgamma(1000, 1, 1), rgamma(1000, 2, 1), rgamma(1000, 1, 2)))
}
