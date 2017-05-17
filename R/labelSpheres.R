labelSpheres <- function(coords, labels, naive=FALSE)
# Spreads the labels to all hyperspheres, based on the closest labelled hypersphere 
# 
# written by Aaron Lun 
# created 15 May 2017
{
    stopifnot(identical(nrow(coords), length(labels)))
    has.label <- labels!=""
    if (!any(has.label) | all(has.label)) { 
        return(labels) 
    }

    fresh.labels <- labels[has.label]
    ulabels <- unique(fresh.labels)
    if (length(ulabels)==1) { 
        labels[] <- ulabels
        return(labels)
    }

    # Preparing data for counting.
    coords <- as.matrix(coords)    
    labelled <- coords[has.label,,drop=FALSE]
    npts <- nrow(labelled)
    if (!naive) { 
        converted <- .reorganize_cells(exprs=labelled, sample.id=rep(1L, npts), cell.id=seq_len(npts))
        label.coords <- converted$exprs
        metadata <- converted$metadata
        hyper.ids <- converted$cell.id
    } else {
        label.coords <- t(labelled)
        hyper.ids <- seq_len(npts)
        metadata <- list()
    }
    cluster.centers <- metadata$cluster.centers
    cluster.info <- metadata$cluster.info

    # Finding closest labelled hypersphere for each other hypersphere.
    closest <- .Call(cxx_find_knn, label.coords, cluster.centers, cluster.info, 1L, -2L, t(coords)) 
    if (is.character(closest)) {
        stop(closest)
    }
    new.labels <- fresh.labels[hyper.ids][closest+1L]
    return(new.labels)
}
