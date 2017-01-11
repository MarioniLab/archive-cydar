\name{dnaGate}
\alias{dnaGate}

\title{Gate events based on DNA channels}
\description{
Construct a gate to remove debris and doublets, based on the two DNA (iridium) channels used in most mass cytometry experiments. 
}

\usage{
dnaGate(x, name1, name2, tol=0.5, nmads=3, type=c("both", "lower"), 
    shoulder=FALSE, rank=1, ...)
}

\arguments{
\item{x}{A flowFrame object like that constructed by \code{\link{poolCells}}.}
\item{name1, name2}{Strings containing the names of the two DNA channels.}
\item{tol}{A numeric scalar quantifying the maximum distance from the equality line.}
\item{nmads}{A numeric scalar specifying the number of median absolute deviations (MADs) beyond which an event can be considered an outlier.}
\item{type}{A string specifying the type of gating to be performed.}
\item{shoulder}{A logical scalar indicating whether the function should attempt to detect shoulders.}
\item{rank}{An integer scalar specifying the peak corresponding to singlets. 
By default, the largest mode is treated as the singlet peak.}
\item{...}{Additional arguments to pass to \code{\link{density}}, to fine-tune identification of local minima.}
}

\details{
For each DNA channel, the \code{rank}th-largest local mode is identified and is assumed to correspond to singlets.
Local minima of density that neighbour the chosen mode are identified. 
To remove debris, the lower bound is set to the largest local minima that is smaller than the chosen mode.
To remove doublets, the upper bound is set to the smallest local minima that is larger than the chosen mode.

We also consider an alternative lower bound at \code{nmads} MADs below the chosen mode.
(Here, the MAD is computed using only values below the mode, to avoid potential inflation due to a doublet peak.)
If this alternative is larger than the largest local minima below the mode, it is used as the lower bound instead.
This avoids using a poor lower bound when there are no obvious minima in the distribution.
Similarly, an alternative upper bound is defined at \code{nmads} MADs above the median, and is used if it is smaller than the smallest local minima above the mode. 

For some data sets, there may not be any clear bimodality in the intensity distribution, e.g., if the mean shift is dominated by noise.
If \code{shoulder=TRUE}, the function will attempt to identify the doublet peak as a \dQuote{shoulder} off the singlet peak.
Alternatively, if there is no evidence for separate singlet/doublet peaks, it may not be feasible (or desirable) to try to distinguish them.
In such cases, users can set \code{type="lower"}, whereby the upper bound is set to an arbitrarily large value and effectively ignored during gating.

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
\code{\link{poolCells}},
\code{\link{density}}
}

\examples{
set.seed(200)

### Mocking up some data with clear bimodality: ###
library(flowCore)
singlets <- rnorm(20000, 2, 0.2)
dna1 <- matrix(rnorm(40000, singlets, 0.1), ncol=2)
doublets <- rnorm(10000, 3, 0.2)
dna2 <- matrix(rnorm(20000, doublets, 0.1), ncol=2)
dna.int <- rbind(dna1, dna2)
colnames(dna.int) <- c("Ir191", "Ir193")
ff <- flowFrame(dna.int)

### Defining the gate: ###
dgate <- dnaGate(ff, "Ir191", "Ir193")
smoothScatter(dna.int[,1], dna.int[,2])
polygon(dgate@boundaries[,1], dgate@boundaries[,2], border="red")

### Mocking up some data with no obvious bimodality: ###
singlets <- rnorm(20000, 2, 0.2)
dna1 <- matrix(rnorm(40000, singlets, 0.1), ncol=2)
doublets <- rnorm(10000, 2.5, 0.2) # <- less separation between modes
dna2 <- matrix(rnorm(20000, doublets, 0.1), ncol=2)
dna.int <- rbind(dna1, dna2)
colnames(dna.int) <- c("Ir191", "Ir193")
ff <- flowFrame(dna.int)

### Defining the gate: ###
dgate <- dnaGate(ff, "Ir191", "Ir193", shoulder=TRUE)
smoothScatter(dna.int[,1], dna.int[,2])
polygon(dgate@boundaries[,1], dgate@boundaries[,2], border="red")
}


