\name{dnaGate}
\alias{dnaGate}

\title{Gate events based on DNA channels}
\description{
Construct a gate to remove debris and doublets, based on the two DNA (iridium) channels used in most mass cytometry experiments. 
}

\usage{
dnaGate(x, name1, name2, tol=0.5, nmads=3, type=c("both", "lower"))
}

\arguments{
\item{x}{A flowFrame object like that constructed by \code{\link{poolCells}}.}
\item{name1, name2}{Strings containing the names of the two DNA channels.}
\item{tol}{A numeric scalar quantifying the maximum distance from the equality line.}
\item{nmads}{A numeric scalar specifying the number of median absolute deviations (MADs) beyond which an event can be considered an outlier.}
\item{type}{A string specifying the type of gating to be performed.}
}

\details{
For each DNA channel, the global mode is identified and is assumed to correspond to singlets.
Local minima of density neighbouring the global mode are then identified. 
To remove debris, the lower bound is set to the largest local minima that is smaller than the global mode.
To remove doublets, the upper bound is set to the smallest local minima that is larger than the mode.

We also consider an alternative lower bound at \code{nmads} MADs below the median of the intensity distribution.
If this is larger than the largest local minima, it is used as the lower bound instead.
This avoids using a poor lower bound when there are no obvious minima in the distribution.
Similarly, an alternative upper bound is defined at \code{nmads} MADs above the median, and is used if it is smaller than the smallest local minima. 

If \code{type="lower"}, the upper bound is set to an arbitrarily large value and effectively ignored during gating.
This should be used in cases where there is no clear bimodality in the intensity distribution.
In such cases, it may not be feasible to distinguish between singlets and doublets.

To simultaneously gate on both DNA channels, we fit a line to the paired intensities for all events, i.e., the \dQuote{equality line}.
Two perpendicular lines passing through the paired lower/upper bounds are constructed.
Two parallel lines that are \code{tol} away from the equality line are also defined. 
The box defined by these four lines is used to construct a polygonGate object, within which all events are retained.

The value of \code{tol} represents the maximum Euclidean distance of any event from the equality line in the two-dimensional space.
Any event more the \code{tol} from the line is removed as the two iridium isotopes have not been evenly captured.
This may be indicative of a problem with the TOF detector for this event.
}

\value{
A polygonGate object, defined to retain singlet events.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{polygonGate}},
\code{\link{poolCells}}
}

\examples{
library(flowCore)
dna1 <- matrix(rnorm(40000, 2, 0.2), ncol=2)
dna2 <- matrix(rnorm(20000, 3, 0.2), ncol=2)
dna.int <- rbind(dna1, dna2)
colnames(dna.int) <- c("Ir191", "Ir193")
ff <- flowFrame(dna.int)

dgate <- dnaGate(ff, "Ir191", "Ir193")
smoothScatter(dna.int[,1], dna.int[,2])
polygon(dgate@boundaries[,1], dgate@boundaries[,2], border="red")

}


