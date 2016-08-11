intensityRanges <- function(cell.data, p=0.01)
# Computes the log-fold change in intensity.
#
# written by Aaron Lun
# created 22 April 2016
# last modified 16 June 2016
{
    .check_cell_data(cell.data)
    marker.names <- attributes(cell.data)$markers
    all.ranges <- list()
    for (m in seq_along(marker.names)) {
        all.ranges[[m]] <- quantile(cell.data[m,], p=c(p, 1-p))
    }
    output <- do.call(cbind, all.ranges)
    colnames(output) <- marker.names
    rownames(output) <- c("min", "max")
    return(output)
}
