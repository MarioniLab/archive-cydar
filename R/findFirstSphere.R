findFirstSphere <- function(coords, threshold=1, naive=FALSE)
# Returns a logical vector indicating which hyperspheres are redundant
# within the specified distance threshold.
#
# written by Aaron Lun
# created 31 October 2016
{
    colnames(coords) <- seq_len(ncol(coords)) # dummy colnames to keep it happy.
    converted <- prepareCellData(list(X=coords))
    if (naive) { 
        cluster.centers <- cluster.info <- NULL
    } else {
        cluster.centers <- attributes(converted)$cluster.centers
        cluster.info <- attributes(converted)$cluster.info
    }

    threshold <- as.double(threshold)
    nred <- .Call(cxx_drop_redundant, converted, cluster.centers, cluster.info, threshold)
    if (is.character(nred)) { stop(nred) }

    return(data.frame(keep=nred[[1]], closest=nred[[2]], distance=nred[[3]]))
}

