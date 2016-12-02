meanIntensities <- function(x, markers) 
# Computes the mean intensity for each sample, in each hypersphere,
# for each marker. The idea is to fit this to a linear model.
#
# written by Aaron Lun
# created 2 December 2016
{
    .check_cell_data(x, check.clusters=FALSE)
    samples <- colnames(x)
    used <- .chosen_markers(markers, markernames(x))
    sample.id <- cellData(x)$sample.id - 1L # Get to zero indexing.

    ci <- cellIntensities(x)
    out <- .Call(cxx_compute_mean_int, ci, length(samples), sample.id, cellAssignments(x), used)
    if (is.character(out)) stop(out)
    
    used.markers <- markernames(x)[used]
    old.names <- assayNames(x)
    for (j in seq_along(used.markers)) { 
        assay(x, length(old.names)+j) <- out[[j]]
    }
    assayNames(x) <- c(old.names, paste0("mean.", used.markers))
    return(x)
}

