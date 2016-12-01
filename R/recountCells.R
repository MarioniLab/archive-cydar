recountCells <- function(x, markers, tol=0.5, filter=10L)
# This function recounts cells in each group in 'x' into a new set of hyperspheres
# in the specified space of 'markers'.
#
# written by Aaron Lun
# created 30 November 2016
# last modified 1 December 2016
{
    all.assign <- cellAssignments(x)
    used <- .chosen_markers(markers, markernames(x))
    ci <- cellIntensities(x)
    distance <- tol*sqrt(sum(used))
    filter <- as.integer(filter)
  
    # Collating values. 
    out <- .Call(cxx_recount_cells, ci, used, distance, rowData(x)$center.cell-1L, all.assign, filter)
    if (is.character(out)) stop(out)
    combined.ass <- unlist(out[[1]], recursive=FALSE)
    combined.centers <- out[[2]]
    combined.groups <- rep(seq_len(nrow(x)), lengths(combined.centers))
    combined.centers <- unlist(combined.centers)

    # Computing assorted statistics.
    nsamples <- ncol(x)
    sample.id <- cellData(x)$sample.id - 1L
    out.stats <- .Call(cxx_compute_hyperstats, ci, nsamples, sample.id, combined.ass, used) 
    if (is.character(out.stats)) stop(out.stats)

    combined.counts <- out.stats[[1]]
    colnames(combined.counts) <- sampleNames(x)
    combined.coords <- out.stats[[2]]
    colnames(combined.coords) <- markernames(x)

    # Returning a cyData object.
    new.markerdata <- markerData(x)
    new.markerdata$reused <- used
    new.metadata <- metadata(x)
    new.metadata$retol <- tol
    cyData(assays=list(counts=combined.counts),
           intensities=combined.coords,
           cellAssignments=combined.ass,
           cellIntensities=cellIntensities(x),
           markerData=new.markerdata,
           cellData=cellData(x),
           rowData=DataFrame(center.cell=combined.centers, group=combined.groups),
           colData=colData(x),
           metadata=new.metadata
           ) 
}

