neighborDistances <- function(cell.data, neighbors=50, downsample=50, as.tol=TRUE, naive=FALSE) 
# Calculates the 'tol' required to capture a certain number of neighbors.
#
# written by Aaron Lun
# created 7 July 2016
# last modified 11 August 2016
{
    .check_cell_data(cell.data)
    markers <- attributes(cell.data)$markers
    if (naive) {
        cluster.centers <- cluster.info <- NULL
    } else {
        cluster.centers <- attributes(cell.data)$cluster.centers
        cluster.info <- attributes(cell.data)$cluster.info
    }

    # Checking parameters.
    neighbors <- as.integer(neighbors)
    if (neighbors <= 0L) { stop("'neighbors' should be a positive integer") }
    downsample <- as.integer(downsample)
    if (downsample <= 0L) { stop("'downsample' should be a positive integer") }

    # Computing distances.
    distances <- .Call(cxx_get_nndist, cell.data, cluster.centers, cluster.info, neighbors, downsample) 
    if (is.character(distances)) {
        stop(distances)
    }

    # Converting to tolerance values, if so desired.
    distances <- t(distances)
    if (as.tol) {
        distances <- distances/sqrt(length(markers))
    }
    return(distances)
}
