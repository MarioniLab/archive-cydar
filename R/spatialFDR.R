spatialFDR <- function(coords, pvalues, neighbors=50, bandwidth=NULL, naive=FALSE) 
# This controls the spatial FDR across a set of plot coordinates.
# Each point is weighted by the reciprocal of its density, based on the specified 'radius'.
# A frequency-weighted version of the BH method is then applied to the p-values.
#
# written by Aaron Lun
# created 23 May 2016
# last modified 15 May 2017
{
    if (length(pvalues)!=nrow(coords)) { stop("coords 'nrow' and p-value vector length are not the same") }
    colnames(coords) <- seq_len(ncol(coords)) # dummy colnames to keep it happy.
    coords <- .find_valid_markers(coords, return.matrix=TRUE) # dropping NA markers.

    # Discarding NA pvalues.
    haspval <- !is.na(pvalues)
    if (!all(haspval)) {
        coords <- coords[haspval,,drop=FALSE]
        pvalues <- pvalues[haspval]
    }

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

    if (is.null(bandwidth)) { 
        neighbors <- as.integer(neighbors)
        if (neighbors==0L) { 
            bandwidth <- 0 
        } else if (neighbors < 0L) { 
            stop("'neighbors' must be a non-negative integer") 
        } else { 
            # Figuring out the bandwidth for KDE, as the median of distances to the n-th neighbour.
            allbands <- .Call(cxx_find_knn, new.coords, cluster.centers, cluster.info, neighbors, -1L, NULL)
            if (is.character(allbands)) { stop(allbands) }
            bandwidth <- median(allbands)
        }
    } else {
        bandwidth <- as.double(bandwidth)
    }
    if (bandwidth <= 0) {
        warning("setting a non-positive bandwidth to a small offset")
        bandwidth <- 1e-8 
    }

    # Computing densities with a tricube kernel.
    densities <- .Call(cxx_compute_density, new.coords, cluster.centers, cluster.info, bandwidth)
    if (is.character(densities)) { stop(densities) }
    w <- 1/densities
    w[hyper.ids] <- w # Getting back to original ordering.

    # Computing a density-weighted q-value.
    o <- order(pvalues)
    pvalues <- pvalues[o]
    w <- w[o]

    adjp <- numeric(length(o))
    adjp[o] <- rev(cummin(rev(sum(w)*pvalues/cumsum(w))))
    adjp <- pmin(adjp, 1)

    if (!all(haspval)) {
        refp <- rep(NA_real_, length(haspval))
        refp[haspval] <- adjp
        adjp <- refp
    }
    return(adjp)
}

.find_valid_markers <- function(coords, return.matrix=FALSE) {
    NC <- ncol(coords)
    to.use <- logical(NC)
    for (i in seq_len(NC)) {
        failed <- is.na(coords[,i])
        if (!any(failed)) {
            to.use[i] <- TRUE
        } else if (!all(failed)) {
            stop("columns should be all NA or with no NA values")
        }
    }
    if (!return.matrix) return(to.use)  
    if (!all(to.use)) coords <- coords[,to.use,drop=FALSE]
    return(coords)
}

