\name{diffIntDistr}
\alias{diffIntDistr}

\title{Compute differences in intensity distributions}

\description{Calculate the difference in the intensity distributions between samples across multiple batches.}

\usage{
diffIntDistr(..., markers=NULL, npts=200) 
}

\arguments{
\item{...}{
One or more lists of intensity matrices or ncdfFlowSet objects, like \code{x} in \code{?\link{prepareCellData}}.
These can be named, otherwise default names will be assigned.
}
\item{markers}{
A character vector specifying which markers should have their differences computed.
All markers that are present in all objects will be used by default.
}
\item{npts}{An integer scalar specifying the number of points for computing differences.}
}

\details{
The difference between two intensity distributions is defined as half of the area between their probability density functions.
For each intensity distribution, the density function is approximated by binning events into a grid of length \code{npts}.
The absolute difference in the proportion of events is computed for each bin; summed across all bins; multiplied by the width of each bin; and halved.
The final value ranges from 0 to 1 and represents the proportion of probability that differs between the distributions.

The differences for each marker are returned in the form of a symmetric matrix.
Each row/column is a sample, named as a combination of the batch and sample names.
Samples are ordered according to the specified order of batches in \code{...} and the order of samples in each batch.
Each entry of the matrix contains the difference between the corresponding samples.

If multiple objects are specified in \code{...}, they are assumed to represent different batches of samples.
The differences between samples in the same or different batches are useful for evaluating the effects of inter-batch normalization, as shown in the example below.
If normalization was successful, differences between batches should be similar to the differences within batches. 
However, this diagnostic is only applicable if the batches are not confounded by biological factors - 
otherwise, if such factors were present, inter-batch differences \emph{should} be larger.
}

\value{
A list is returned with two objects:
\describe{
\item{\code{index}:}{
An integer matrix specifying which entries of each matrix correspond to differences between samples in the same batch.
These are numbered according to the specified order of the batches in \code{...}.
Between-batch differences have values of 0.
}
\item{\code{difference}:}{
A named list of symmetric numeric matrices, with one for each marker.
Each row and column represents a sample and is named accordingly. 
Each entry of the matrix contains the difference between intensity distributions for the two samples corresponding to the row and column.
}
}
}

\seealso{
\code{\link{normalizeBatch}}
}

\author{
Aaron Lun
}

\examples{
example(normalizeBatch)

original <- do.call(diffIntDistr, all.x)
par(mfrow=c(5,4), mar=c(2.1, 2.1, 2.1, 1.1))
for (m in names(original$difference)) { 
    current <- original$difference[[m]] * 100
    is.inter <- original$index==0L
    is.intra <- !is.inter & lower.tri(original$index)
    boxplot(list(Inter=current[is.inter], Intra=current[is.intra]), main=m)
}

par(mfrow=c(5,4), mar=c(2.1, 2.1, 2.1, 1.1))
fixed <- do.call(diffIntDistr, corrected)
for (m in names(fixed$difference)) { 
    current <- fixed$difference[[m]] * 100
    is.inter <- fixed$index==0L
    is.intra <- !is.inter & lower.tri(fixed$index)
    boxplot(list(Inter=current[is.inter], Intra=current[is.intra]), main=m)
}


}
