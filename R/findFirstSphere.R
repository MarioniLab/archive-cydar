findFirstSphere <- function(coords, pvalues, threshold=1, naive=FALSE)
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

    # Getting the original ordering.
    cell.ids <- attributes(converted)$cell.id + 1L

    # Checking for non-redundancy.
    threshold <- as.double(threshold)
    actual_order <- order(pvalues[cell.ids]) - 1L
    out <- .Call(cxx_drop_redundant, actual_order, converted, cluster.centers, cluster.info, threshold)
    if (is.character(out)) { stop(out) }
   
    # Generating output with the original ordering in 'coords' 
    out[cell.ids] <- out
    return(out)
}

