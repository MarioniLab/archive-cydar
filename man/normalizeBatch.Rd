\name{normalizeBatch}
\alias{normalizeBatch}

\title{Normalize intensities across batches}
\description{Perform normalization to correct intensities across batches with at least one common level.}

\usage{
normalizeBatch(batch.x, batch.comp, mode=c("range", "curve"), 
    p=0.01, ref=1, npts=512, extra.init=NULL)
}

\arguments{
\item{batch.x}{A list, where each element is of the same type as \code{x} used in \code{\link{prepareCellData}} (i.e., a ncdfFlowSet or a list of intensity matrices across all samples).}
\item{batch.comp}{
A list of factors (or elements coercible to factors) specifying the composition of each batch, i.e., which samples belong to which groups.
Also can be \code{NULL}, see below.
}
\item{mode}{A string specifying whether range- or curve-based normalization should be performed.}
\item{p}{A numeric scalar between 0 and 0.5, specifying the percentile used to define the range of the distribution for range-based normalization.}
\item{ref}{An integer scalar specifying which batch is to be used as the reference for curve-based normalization.}
\item{npts}{An integer scalar specifying the number of points for the area calculation in curve-based normalization, passed to \code{\link{density}}.}
\item{extra.init}{A list of named numeric vectors for use in \code{\link{optim}} during curve-based normalization, see below for details.}
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
\item If \code{mode="curve"}, an adjustment function is constructed that minimizes the area between probability density functions (PDFs) of the pooled intensity distributions.
The \code{ref} batch is considered to be the reference, and the pooled distributions of each other batches are adjusted towards it using \code{\link{optim}}.
The value of \code{npts} determines the resolution of the area calculation while \code{extra.init} specifies additional starting points to test during optimization.
See below for more details on how curve-based normalization works.
\item The scaling or adjustment function is applied to the intensities of all samples in that batch, yielding corrected intensities for direct comparisons between samples.
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

\section{Details on curve-based normalization}{ 
This approach assumes that the batch effect involves a monotonic transformation of intensities.
In particular, a modified version of the logistic function is considered:

\deqn{
\frac{A}{1 +\exp[B-x\exp(C)]} - \frac{A}{1 + \exp(B)}
}

where \eqn{x} is the intensity, \eqn{y} is the adjustment to the intensity, and the other terms are various constants to be estimated. 
The normalized intensity is simply \code{x+y}.
This function was chosen as it passes through the origin, so that zero intensities are unaffected by normalization;
it can perform a sharp transition from no adjustment to a large adjustment, if only high intensities are distorted between batches;
and it can mimic a straight line, to rescale peak widths between batches.

The terms \eqn{A}, \eqn{B} and \eqn{C} are chosen to minimize the area between the PDFs of the pooled distributons of different batches.
This area is calculated based on density estimates from \code{\link{density}}, with resolution determined by \code{npts}.
Zero intensities are handled separately as they do not form a continuum.
The aim is to adjust the distributions so that they become as similar as possible to the distribution of the reference batch, as specified in \code{ref}.

Minimization is performed using \code{\link{optim}} from a variety of initiation points.
Advanced users can specify additional starting points in \code{extra.init} if the algorithm becomes trapped in some local minima.
This should be a list of named numeric vectors, where each vector contains \code{A}, which controls the direction and magnitude of the correction;
\code{B}, which controls the location of the transition point; and \code{C}, which controls the gradient of the transition.
}

\section{Choosing between normalization methods}{
Curve-based normalization is more powerful than range-based normalization, as the former can eliminate non-linear changes to the intensities whereas the latter cannot.
However, it expects that the shape of the pooled distribution is roughly similar across batches.
Large differences (e.g., a peak present in one batch and absent in another) may lead to incorrect adjustments.

Such differences may be present when batches are confounded with uninteresting biological factors (e.g., individual, mouse of origin) that affect cell abundance. 
In such cases, range-based normalization with \code{mode="range"} is recommended as it is more constrained in how the intensities are adjusted.
This reduces the risk of distorting the intensities, albeit at the cost of \dQuote{under-normalizing} the data.

See also \code{\link{diffIntDist}} for a method to assess the effects of batch normalization.
}

\value{
A list of lists, where each internal list corresponds to a batch and contains intensity matrices corresponding to all samples in that batch.
This matches the format of \code{batch.x}.
}

\author{Aaron Lun}

\seealso{
\code{\link{prepareCellData}},
\code{\link{diffIntDist}}
}

\examples{
### Mocking up some data: ###
nmarkers <- 20
marker.names <- paste0("X", seq_len(nmarkers))
all.x <- list()

for (b in paste0("Batch", 1:3)) { # 3 batches
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
corrected <- normalizeBatch(all.x, batch.comp, mode="curve")
par(mfrow=c(1,2))
plot(ecdf(all.x[[1]][[3]][,1]), col="blue", main="Before")
plot(ecdf(all.x[[2]][[3]][,1]), add=TRUE, col="red")
plot(ecdf(corrected[[1]][[3]][,1]), col="blue", main="After")
plot(ecdf(corrected[[2]][[3]][,1]), add=TRUE, col="red")
}

