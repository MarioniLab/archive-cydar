\name{normalizeBatch}
\alias{normalizeBatch}

\title{Quantile normalize intensities across batches}
\description{Perform quantile normalization to correct intensities across batches with at least one common level.}

\usage{
normalizeBatch(batch.x, batch.comp, mode=c("quantile", "range"), p=0.01)
}

\arguments{
\item{batch.x}{A list, where each element is of the same type as \code{x} used in \code{\link{prepareCellData}} (i.e., a ncdfFlowSet or a list of intensity matrices across all samples).}
\item{batch.comp}{A list of factors (or elements coercible to factors) specifying the composition of each batch, i.e., which samples belong to which groups.}
\item{mode}{A string specifying whether quantile or range-based normalization should be performed.}
\item{p}{A numeric scalar between 0 and 0.5, specifying the percentile used to define the range of the distribution for range-based normalization.}
}

\details{ 
Consider an experiment containing several batches of barcoded samples, in which the barcoding was performed within but not between batches.
This function performs quantile normalization on the intensities for each marker such that they are comparable between samples in different batches.
It relies on the presence of samples from the same group across the different batches to construct the adjustment function for each batch.
This function is then applied to all other samples in each batch, and the entire process is repeated with the next marker.

The process for each marker is as follows. 
Firstly, a quantile function is constructed for each batch, from the weighted pool of intensities from all samples of that batch (see below for details about weighting).
These functions are averaged across batches to obtain a reference quantile function, representing a reference distribution.
For each batch, an adjustment function is constructed that converts the batch-specific quantile function to that of the reference quantile function.
This adjustment is applied to the intensities of all samples in that batch, yielding corrected intensities for direct comparisons between samples.

Weighting is performed to downweight the contribution of larger samples within each batch, as well as to match the composition of different batches.
For example, a batch with more samples in group A and fewer samples in group B would downweight the former and upweight the latter.
Groupings can be specified as batch-specific factors in \code{batch.comp}, with at least one common group required across all batches.
Ideally, all batches would contain samples from all groups, with similar total numbers of cells across batches for each group.
For non-shared groups, note that the adjustment function may not be accurate if the to-be-corrected intensities lie outside the range of values used to construct the function.

If \code{mode="range"}, range-based normalization is performed instead instead of quantile normalization.
For each marker, the range of the reference distribution is computed, using percentiles at \code{p} to avoid distortions due to outliers.
A batch-specific scaling function is defined to equalize the range of the weighted distribution of intensities from each batch to the reference range.
This function is then applied to all samples, yielding normalized intensities.
Compared to \code{mode="quantile"}, range-based normalization avoids assuming that all differences between batches are technical, but cannot handle non-linear technical effects.

To convert the output into a format appropriate for \code{\link{prepareCellData}}, apply \code{\link{unlist}} with \code{recursive=FALSE}.
This will generate a list of intensity matrices for all samples in all batches, rather than a list of list of matrices.
Note that a batch effect should still be included in the design matrix when modelling abundances, as only the intensities are corrected here.
}

\value{
A list of lists, where each internal list corresponds to a batch and contains intensity matrices corresponding to all samples in that batch.
This matches the format of \code{batch.x}.
}

\author{Aaron Lun}

\seealso{
\code{\link{prepareCellData}}
}

\examples{
### Mocking up some data: ###
nmarkers <- 40
marker.names <- paste0("X", seq_len(nmarkers))
all.x <- list()

for (b in 1:3) { # 3 batches
    nsamples <- 10
    sample.names <- paste0("Y", seq_len(nsamples))
    trans.shift <- rnorm(nmarkers, 0, 1)
    trans.grad <- runif(nmarkers, 1, 2)
    x <- list()
    for (i in sample.names) {
        ex <- matrix(rgamma(nmarkers*1000, 2, 2), ncol=nmarkers, nrow=1000)
        ex <- t(t(ex)*trans.grad + trans.shift)
        colnames(ex) <- marker.names
        x[[i]] <- ex
    }   
    all.x[[b]] <- x
}

batch.comp <- list( # Each batch contains different composition/ordering of groups
    factor(rep(1:2, c(3,7))),
    factor(rep(1:2, c(7,3))),
    factor(rep(1:2, 5))
)

### Running the function: ###
corrected <- normalizeBatch(all.x, batch.comp)
par(mfrow=c(1,2))
plot(ecdf(all.x[[1]][[3]][,1]), col="blue", main="Before")
plot(ecdf(all.x[[2]][[3]][,1]), add=TRUE, col="red")
plot(ecdf(corrected[[1]][[3]][,1]), col="blue", main="After")
plot(ecdf(corrected[[2]][[3]][,1]), add=TRUE, col="red")

# Similar effects with range-based normalization.
# (Only in this simple example with scaling and shifting.)
corrected <- normalizeBatch(all.x, batch.comp, mode="range")
par(mfrow=c(1,2))
plot(ecdf(all.x[[1]][[3]][,1]), col="blue", main="Before")
plot(ecdf(all.x[[2]][[3]][,1]), add=TRUE, col="red")
plot(ecdf(corrected[[1]][[3]][,1]), col="blue", main="After")
plot(ecdf(corrected[[2]][[3]][,1]), add=TRUE, col="red")
}

