\name{labelSpheres}
\alias{labelSpheres}

\title{Label unannotated hyperspheres}
\description{Given a set of labels for annotated hyperspheres, propagate labels to the surrounding unannotated hyperspheres.}

\usage{
labelSpheres(coords, labels, naive=FALSE)
}

\arguments{
\item{coords}{A numeric matrix of hypersphere coordinates, containing the median intensity of each marker (column) in each hypersphere (row).}
\item{labels}{A character vector of labels for each hypersphere, set to an empty string for unannotated hyperspheres.}
\item{naive}{A logical scalar specifying whether a naive search should be performed.}
}

\details{
After some hyperspheres have been labelled with \code{\link{interpretSpheres}}, the remainder can be automatically labelled with this function.
Unlabelled hyperspheres are assigned the label of the closest labelled hypersphere.
Obviously, this assumes that enough hyperspheres have been labelled so that the closest hypersphere is of a similar cell type/state.
}

\value{
A character vector containing labels for all hyperspheres.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{interpretSpheres}}
}

\examples{
set.seed(1000)
coords <- matrix(rgamma(10000, 2, 2), nrow=1000)
labels <- character(nrow(coords))
labels[1:4] <- c("B", "CD4T", "CD8T", "Mono")

ref <- labelSpheres(coords, labels)
naive <- labelSpheres(coords, labels, naive=TRUE)

\dontshow{
stopifnot(identical(ref, naive))
}
}

