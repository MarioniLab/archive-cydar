\name{outlierGate}
\alias{outlierGate}

\title{Create an outlier gate}

\description{
Define gating thresholds to remove outlier events for a particular channel.
}

\usage{
outlierGate(x, name, nmads=3, type=c("both", "upper", "lower"))
}

\arguments{
\item{x}{A flowFrame object like that constructed by \code{\link{poolCells}}.}
\item{name}{A string specifying the name of the channel in \code{x} from which intensities are to be extracted.}
\item{nmads}{A numeric scalar specifying the number of median absolute deviations (MADs) beyond which an event can be considered an outlier.}
\item{type}{A string specifying the type of outliers to be removed.}
}

\details{
Outliers are defined as events with intensities that are more than \code{nmads} median absolute deviations from the median of the intensity distribution.
The lower gate threshold is defined as the median minus \code{nmads} MADs, while the upper gate threshold is defined as the median plus \code{nmads} MADs.
If \code{type="upper"}, only large outliers are removed (e.g., dead/alive stains), so the lower threshold is set to \code{-Inf}.
If \code{type="lower"}, only small outliers are removed (e.g., DNA), so the upper threshold is set to \code{Inf}.
}

\value{
A rectangleGate object with lower and upper thresholds defined from \code{x}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{poolCells}},
\code{\link[flowCore]{rectangleGate}}
}

\examples{
example(poolCells)
ogate <- outlierGate(ff, "X1")
ogate

ogate <- outlierGate(ff, "X2", type="upper")
ogate

ogate <- outlierGate(ff, "X3", type="lower")
ogate

sff <- Subset(ff, ogate) # for actual gating.
}
