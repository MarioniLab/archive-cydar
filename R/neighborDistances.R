neighborDistances <- function(x, neighbors=50, downsample=50, as.tol=TRUE, naive=FALSE) 
# Calculates the 'tol' required to capture a certain number of neighbors.
#
# written by Aaron Lun
# created 7 July 2016
# last modified 29 November 2016
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
    nmarkers <- sum(used)

    # Checking parameters.
    neighbors <- as.integer(neighbors)
    if (neighbors <= 0L) { stop("'neighbors' should be a positive integer") }
    downsample <- as.integer(downsample)
    if (downsample <= 0L) { stop("'downsample' should be a positive integer") }

    # Computing distances.
    distances <- .Call(cxx_get_nndist, cellIntensities(x), used, cluster.centers, cluster.info, neighbors, downsample) 
    if (is.character(distances)) {
        stop(distances)
    }

    # Converting to tolerance values, if so desired.
    distances <- t(distances)
    if (as.tol) {
        distances <- distances/sqrt(nmarkers)
    }
    return(distances)
}
