\name{packIndices}
\alias{packIndices}
\alias{unpackIndices}

\title{Pack or unpack indices}
\description{Compress or decompress a vector of indices, usually specifying cell assignments to groups.}

\usage{
packIndices(assignments)
unpackIndices(assignments)
}

\arguments{
\item{assignments}{A list of integer vectors containing indices that specify cell assignments to groups.
Each vector corresponds to a group and should be uncompressed for \code{packIndices}, or already compressed for \code{unpackIndices}.
}
}

\details{
Indices are stored in a compressed format whereby a negative number indicates that all consecutive integers from the preceding index should be used.
For example, a sequence of \code{c(1, 3, -6, 8, -10, 12)} would be unpacked as \code{c(1, 3:6, 8:10, 12)}.
This saves a lot of memory for storing cell assignments in hyperspheres, where there are likely to be many consecutive indices.

The \code{\link{countCells}} function will automatically fill in the \code{cellAssignments} slot of the output object with compressed index vectors.
The \code{unpackIndices} function can be used to construct the full vector from some or all of these vectors, for manual use elsewhere.
Conversely, users manually constructing cyData objects can use \code{packIndices} to compress the vectors and save space.
}

\value{
A list of integer vectors containing compressed or uncompressed index vectors, for \code{packIndices} and \code{unpackIndices} respectively.
}

\seealso{
\code{\link{countCells}}
}

\author{Aaron Lun}

\examples{
a <- c(1L, 3:6, 8:10, 12L)
packIndices(list(a))

b <- c(1L, 3L, -6L, 8L, -10L, 12L)
unpackIndices(list(b))

}
