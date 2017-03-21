intensityRanges <- function(x, p=0.01)
# Computes the log-fold change in intensity.
#
# written by Aaron Lun
# created 22 April 2016
# last modified 21 March 2017
{
    .check_cell_data(x, check.clusters=FALSE)
    marker.names <- rownames(markerData(x))
    ci <- cellIntensities(x)
    all.ranges <- vector("list", length(marker.names))
    for (m in seq_along(marker.names)) {
        all.ranges[[m]] <- quantile(ci[m,], p=c(p, 1-p))
    }
    output <- do.call(cbind, all.ranges)
    colnames(output) <- marker.names
    rownames(output) <- c("min", "max")
    return(output)
}
