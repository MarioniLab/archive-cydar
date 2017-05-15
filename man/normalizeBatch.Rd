\name{normalizeBatch}
\alias{normalizeBatch}

\title{Normalize intensities across batches}
\description{Perform normalization to correct intensities across batches with at least one common level.}

\usage{
normalizeBatch(batch.x, batch.comp, mode="range", p=0.01)
}

\arguments{
\item{batch.x}{A list, where each element is of the same type as \code{x} used in \code{\link{prepareCellData}} (i.e., a ncdfFlowSet or a list of intensity matrices across all samples).}
\item{batch.comp}{
A list of factors (or elements coercible to factors) specifying the composition of each batch, i.e., which samples belong to which groups.
Also can be \code{NULL}, see below.
}
\item{mode}{
A string or character vector of length equal to the number of markers, specifying whether range-based or warping normalization should be performed for each marker.
This can take values of \code{"range"}, \code{"warp"} or \code{"none"} (in which case no normalization is performed).
}
\item{p}{A numeric scalar between 0 and 0.5, specifying the percentile used to define the range of the distribution for range-based normalization.}
}

\details{ 
Consider an experiment containing several batches of barcoded samples, in which the barcoding was performed within but not between batches.
This function normalizes the intensities for each marker such that they are comparable between samples in different batches.
The process for each marker is as follows: 
\enumerate{
\item Weighting is performed to downweight the contribution of larger samples within each batch, as well as to match the composition of samples across different batches.
The composition of each batch can be specified by \code{batch.comp}, see below for more details.
The weighted intensities for each batch represents the pooled distribution of intensities from all samples in that batch.
\item If \code{mode="range"}, a quantile function is constructed for the pooled distribution of each batch.
These functions are averaged across batches to obtain a reference quantile function, representing a reference distribution.
The range of the reference distribution is computed at percentiles \code{p} and \code{1-p} (to avoid distortions due to outliers).
A batch-specific scaling function is defined to equalize the range of the weighted distribution of intensities from each batch to the reference range.
\item If \code{mode="warp"}, weighted sampling from each pooled distribution is performed to generate a pseudo-sample for each batch.
This is used to construct a flowSet for use in warping normalization - see \code{?\link{normalization}} and \code{?\link[flowStats]{warpSet}} for details.
A warping function is computed for each batch that adjust the intensity distribution to be more similar to a reference (i.e., the first batch).
\item The scaling or warping function is applied to the intensities of all samples in that batch, yielding corrected intensities for direct comparisons between samples.
}

Groupings can be specified as batch-specific factors in \code{batch.comp}, with at least one common group required across all batches.
If the composition of each batch is the same, \code{batch.comp} can be set to \code{NULL} rather than being manually specified.
This composition is used to weight the contribution of each sample to the reference distribution.
For example, a batch with more samples in group A and fewer samples in group B would get lower weights assigned to the former and larger weights to the latter.

Construction of the adjustment function relies on the presence of samples from the same group across the different batches.
Ideally, all batches would contain samples from all groups, with similar total numbers of cells across batches for each group.
The adjustment function will still be applied to intensities for samples from non-shared groups that do not contribute to the reference distribution.
However, note that the adjustment may not be accurate if the to-be-corrected intensities lie outside the range of values used to construct the function.

To convert the output into a format appropriate for \code{\link{prepareCellData}}, apply \code{\link{unlist}} with \code{recursive=FALSE}.
This will generate a list of intensity matrices for all samples in all batches, rather than a list of list of matrices.
Note that a batch effect should still be included in the design matrix when modelling abundances, as only the intensities are corrected here.
}

\section{Choosing between normalization methods}{
Warping normalization can be more powerful than range-based normalization, as the former can eliminate non-linear changes to the intensities whereas the latter cannot.
However, it requires that landmarks in the intensity distribution (i.e., peaks) be easily identifiable and consistent across batches. 
Large differences (e.g., a peak present in one batch and absent in another) may lead to incorrect adjustments.

Such differences may be present when batches are confounded with uninteresting biological factors (e.g., individual, mouse of origin) that affect cell abundance. 
In such cases, range-based normalization with \code{mode="range"} is recommended as it is more constrained in how the intensities are adjusted.
This reduces the risk of distorting the intensities, albeit at the cost of \dQuote{under-normalizing} the data.

It is advisable to inspect the intensity distributions before and after normalization, to ensure that the methods have behaved appropriately.
This can be done by constructing histograms for each marker with \code{\link{multiIntHist}}.
See also \code{\link{diffIntDistr}} for a quantitative measure of similarity between distributions.
}

\value{
A list of lists, where each internal list corresponds to a batch and contains intensity matrices corresponding to all samples in that batch.
This matches the format of \code{batch.x}.
}

\author{Aaron Lun}

\seealso{
\code{\link{prepareCellData}},
\code{\link{diffIntDistr}},
\code{\link{multiIntHist}},
\code{\link[flowCore]{normalization}},
\code{\link[flowStats]{warpSet}} 
}

\examples{
### Mocking up some data: ###
nmarkers <- 10
marker.names <- paste0("X", seq_len(nmarkers))
all.x <- list()

for (b in paste0("Batch", 1:3)) { # 3 batches
    nsamples <- 10
    sample.names <- paste0("Y", seq_len(nsamples))
    trans.shift <- runif(nmarkers, 0, 1)
    trans.grad <- runif(nmarkers, 1, 2)
    x <- list()
    for (i in sample.names) {
        ex <- matrix(rgamma(nmarkers*1000, 2, 2), nrow=nmarkers)
        ex <- t(ex*trans.grad + trans.shift)
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
corrected <- normalizeBatch(all.x, batch.comp, mode="range")
par(mfrow=c(1,2))
plot(ecdf(all.x[[1]][[3]][,1]), col="blue", main="Before")
plot(ecdf(all.x[[2]][[3]][,1]), add=TRUE, col="red")
plot(ecdf(corrected[[1]][[3]][,1]), col="blue", main="After")
plot(ecdf(corrected[[2]][[3]][,1]), add=TRUE, col="red")

# Similar effects with warping normalization. 
if (.Platform$OS.type!="windows") {
    wcorrected <- normalizeBatch(all.x, batch.comp, mode="warp")
}
}

