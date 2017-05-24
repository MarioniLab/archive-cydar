recountCells <- function(x, markers, tol=0.5)
# This function recounts cells in each group in 'x' into a new set of hyperspheres
# in the specified space of 'markers'.
#
# written by Aaron Lun
# created 30 November 2016
# last modified 24 May 2017
{
    old.used <- .chosen_markers(markerData(x)$used, markernames(x))
    old.distance <- metadata(x)$tol*sqrt(sum(old.used)) # distance used to construct 'x'.
    used <- old.used | .chosen_markers(markers, markernames(x))
    sq.total.nused <- sqrt(sum(used))
    distance <- tol*sq.total.nused # combined distance.
    eff.tol <- old.distance/sq.total.nused
    if (distance >  old.distance) {
        warning("effective tol is ", round(eff.tol, 3), ", not the specified ", round(tol, 3))
        distance <- old.distance # because we can't count things that were missed in original counting of 'x'.
    }
   
    # Collating values. 
    all.assign <- cellAssignments(x)
    ci <- .get_used_intensities(x, used)
    out <- .Call(cxx_recount_cells, ci, distance, rowData(x)$center.cell-1L, all.assign)
    if (is.character(out)) stop(out)
    combined.ass <- out[[1]]

    # Computing assorted statistics.
    nsamples <- ncol(x)
    sample.id <- cellData(x)$sample.id - 1L
    out.stats <- .Call(cxx_compute_hyperstats, ci, nsamples, sample.id, combined.ass) 
    if (is.character(out.stats)) stop(out.stats)

    combined.counts <- out.stats[[1]]
    colnames(combined.counts) <- sampleNames(x)
    combined.coords <- matrix(NA_real_, nrow(x), nmarkers(x))
    combined.coords[,used] <- out.stats[[2]]
    colnames(combined.coords) <- markernames(x)

    # Returning a CyData object.
    new.markerdata <- markerData(x)
    new.markerdata$used <- used 
    new.metadata <- metadata(x)
    new.metadata$cluster.centers <- new.metadata$cluster.info <- NULL # as the speed-ups aren't applicable with new markers.
    new.metadata$tol <- eff.tol
    CyData(assays=list(counts=combined.counts),
           intensities=combined.coords,
           cellAssignments=combined.ass,
           cellIntensities=cellIntensities(x),
           markerData=new.markerdata,
           cellData=cellData(x),
           rowData=rowData(x),
           colData=colData(x),
           metadata=new.metadata
           ) 
}

