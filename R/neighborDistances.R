neighborDistances <- function(x, neighbors=50, downsample=50, as.tol=TRUE, naive=FALSE) 
# Calculates the 'tol' required to capture a certain number of neighbors.
#
# written by Aaron Lun
# created 7 July 2016
# last modified 24 May 2017
{
    .check_cell_data(x, check.clusters=!naive)
    if (naive) {
        cluster.centers <- cluster.info <- NULL
    } else {
        cluster.centers <- metadata(x)$cluster.centers
        cluster.info <- metadata(x)$cluster.info
    }

    # Calculating the number of used markers.
    used <- .chosen_markers(markerData(x)$used, markernames(x))
    nused <- sum(used)

    # Checking parameters.
    neighbors <- as.integer(neighbors)
    if (neighbors <= 0L) { stop("'neighbors' should be a positive integer") }
    downsample <- as.integer(downsample)
    if (downsample <= 0L) { stop("'downsample' should be a positive integer") }

    # Computing distances.
    ci <- .get_used_intensities(x, used)
    distances <- .Call(cxx_get_nndist, ci, cluster.centers, cluster.info, neighbors, downsample) 
    if (is.character(distances)) {
        stop(distances)
    }

    # Converting to tolerance values, if so desired.
    distances <- t(distances)
    if (as.tol) {
        distances <- distances/sqrt(nused)
    }
    return(distances)
}
