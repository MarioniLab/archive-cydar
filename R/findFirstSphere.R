findFirstSphere <- function(coords, pvalues, threshold=1, block=NULL, naive=FALSE)
# Returns a logical vector indicating which hyperspheres are redundant
# within the specified distance threshold.
#
# written by Aaron Lun
# created 31 October 2016
# last modified 2 December 2016
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
    coords <- .find_valid_markers(coords, return.matrix=TRUE) # dropping NA markers.

    # Preparing data for counting.
    coords <- as.matrix(coords)
    npts <- nrow(coords)
    if (!naive) { 
        converted <- .reorganize_cells(exprs=coords, sample.id=rep(1L, npts), cell.id=seq_len(npts))
        new.coords <- converted$exprs
        metadata <- converted$metadata
        hyper.ids <- converted$cell.id
    } else {
        new.coords <- t(coords)
        hyper.ids <- seq_len(npts)
        metadata <- list()
    }
    cluster.centers <- metadata$cluster.centers
    cluster.info <- metadata$cluster.info

    # Checking for non-redundancy.
    threshold <- as.double(threshold)
    actual_order <- order(pvalues[hyper.ids]) - 1L
    out <- .Call(cxx_drop_redundant, actual_order, new.coords, cluster.centers, cluster.info, threshold)
    if (is.character(out)) { stop(out) }
   
    # Generating output with the original ordering in 'coords' 
    out[hyper.ids] <- out
    return(out)
}

