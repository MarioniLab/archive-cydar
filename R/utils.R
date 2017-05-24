.check_cell_data <- function(x, check.clusters=TRUE) 
# Checks incoming cell data, that it was properly processed by prepareCellData.
{
    sample.id <- cellData(x)$sample.id
    stopifnot(all(sample.id > 0L & sample.id <= ncol(x)))

    central.id <- rowData(x)$center.cell
    if (!is.null(central.id)) { 
        stopifnot(all(central.id > 0L & central.id <= ncells(x)))
    }

    if (check.clusters) {
        cluster.centers <- metadata(x)$cluster.centers        
        if (is.null(cluster.centers)) {
            stop("'cluster.centers' must be defined for non-naive counting")
        }
        stopifnot(nrow(cluster.centers)==sum(markerData(x)$used))
        
        cluster.info <- metadata(x)$cluster.info
        if (is.null(cluster.info)) {
            stop("'cluster.info' must be defined for non-naive counting")
        }
        
        stopifnot(ncol(cluster.centers)==length(cluster.info))
        for (clust in cluster.info) {
            stopifnot(!is.unsorted(clust[[2]]))
            stopifnot(clust[[1]] >= 0L & clust[[1]]+length(clust[[2]]) <= ncells(x))
        }
    }

    invisible(NULL)
}

.chosen_markers <- function(chosen.markers, all.markers) 
# Identifying the markers that were chosen for use.
{
    if (!is.null(chosen.markers)) {
        used <- logical(length(all.markers))
        if (is.character(chosen.markers)) {
            chosen.markers <- match(chosen.markers, all.markers)
            if (any(is.na(chosen.markers))) {
                stop("specified 'markers' not in available set of markers")
            }
        }
        used[chosen.markers] <- TRUE
    } else {
        used <- rep(TRUE, length(all.markers))
    }
    return(used)
}

.get_used_intensities <- function(x, used) 
# Getting the intensities that were used from the CyData object.
{
    ci <- cellIntensities(x)
    if (!all(used)) {
        ci <- ci[used,,drop=FALSE]
    }
    return(ci)
}
    
