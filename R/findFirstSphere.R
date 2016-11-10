findFirstSphere <- function(coords, pvalues, threshold=1, block=NULL, naive=FALSE)
# Returns a logical vector indicating which hyperspheres are redundant
# within the specified distance threshold.
#
# written by Aaron Lun
# created 31 October 2016
# last modified 10 November 2016
{
    if (length(pvalues)!=nrow(coords)) {
        stop("length of 'pvalues' must equal number of rows in 'coords'")
    }
    if (!is.null(block)) {
        # Identifying unique elements within each block.
        if (length(block)!=nrow(coords)) {
            stop("length of 'block' must equal number of rows in 'coords'")
        }
        by.block <- split(seq_along(block), block)
        total.out <- logical(length(block))
        for (b in by.block) {
            total.out[b] <- Recall(coords[b,,drop=FALSE], pvalues[b], threshold=threshold, block=NULL, naive=naive)
        }
        return(total.out)
    }

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

