medIntensities <- function(x, markers) 
# Computes the median intensity for each sample, in each hypersphere,
# for each marker. The idea is to fit this to a linear model.
#
# written by Aaron Lun
# created 2 December 2016
# last modified 3 December 2016
{
    .check_cell_data(x, check.clusters=FALSE)
    samples <- colnames(x)
    used <- .chosen_markers(markers, markernames(x))
    sample.id <- cellData(x)$sample.id - 1L # Get to zero indexing.

    ci <- cellIntensities(x)
    out <- .Call(cxx_compute_median_int, ci, length(samples), sample.id, cellAssignments(x), used)
    if (is.character(out)) stop(out)
    
    used.markers <- markernames(x)[used]
    old.names <- assayNames(x)
    for (j in seq_along(used.markers)) { 
        assay(x, length(old.names)+j) <- out[[j]]
    }
    assayNames(x) <- c(old.names, paste0("med.", used.markers))
    return(x)
}

