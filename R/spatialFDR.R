spatialFDR <- function(coords, pvalues, neighbors=50, bandwidth=NULL, naive=FALSE) 
# This controls the spatial FDR across a set of plot coordinates.
# Each point is weighted by the reciprocal of its density, based on the specified 'radius'.
# A frequency-weighted version of the BH method is then applied to the p-values.
#
# written by Aaron Lun
# created 23 May 2016
# last modified 11 August 2016
{
    if (length(pvalues)!=nrow(coords)) { stop("coords 'nrow' and p-value vector length are not the same") }
    if (naive) {
        new.coords <- t(coords)
        cluster.centers <- cluster.info <- NULL
        hyper.ids <- seq_len(nrow(coords))
    } else {
        colnames(coords) <- seq_len(ncol(coords)) # dummy colnames to keep it happy.
        converted <- prepareCellData(list(X=coords))
        new.coords <- cellIntensities(converted)
        cluster.centers <- metadata(converted)$cluster.centers
        cluster.info <- metadata(converted)$cluster.info
        hyper.ids <- cellData(converted)$cell.id
    }

    if (is.null(bandwidth)) { 
        neighbors <- as.integer(neighbors)
        if (neighbors <= 0L) { stop("'neighbors' must be a positive integer") }

        # Figuring out the bandwidth for KDE, as the median of distances to the 50th neighbour.
        allbands <- .Call(cxx_get_knn_distance, new.coords, cluster.centers, cluster.info, neighbors)
        if (is.character(allbands)) { stop(allbands) }
        bandwidth <- median(allbands)
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
    pmin(adjp, 1)
}

