recountCells <- function(x, markers, tol=0.5, downsample=NULL, filter=10, naive=FALSE)
# This function recounts cells in each group in 'x' into a new set of hyperspheres
# in the specified space of 'markers'.
#
# written by Aaron Lun
# created 30 November 2016
{
    all.assign <- cellAssignments(x)
    used <- .chosen_markers(markers, markernames(x))
    my.markers <- cellIntensities(x)[used,,drop=FALSE]

    my.samples <- sampleNames(x)
    nsamples <- length(my.samples)
    sample.ids <- cellData(x)$sample.id

    # Picking the downsampling frequency.
    ncells <- rowSums(assay(x))
    if (is.null(downsample)) { 
        downsampling <- ceiling(ncells*10/(ncells + 50)) # converges to downsample=10 past 50 cells.
    } else {
        downsample <- as.integer(downsample)
        downsampling <- rep(downsample, length.out=nrow(x))
    }
    
    # Other constants required.
    nmarkers <- nrow(my.markers)
    distance <- tol*sqrt(nmarkers)
    use.markers <- rep(TRUE, nmarkers)

    # Unpacking and recounting.
    combined <- list()
    for (r in seq_along(all.assign)) {
        idx <- unpackIndices(all.assign[r])[[1]]
        cur.data <- my.markers[,idx,drop=FALSE]
        sample.id <- sample.ids[idx]
        cell.id <- seq_along(idx)

        # Preparing the data.
        if (length(idx) <= 20 || naive) {
            metadata <- list()
        } else {
            reorg <- .reorganize_cells(t(cur.data), sample.id, cell.id)   
            sample.id <- reorg$sample.id
            cell.id <- reorg$cell.id
            metadata <- reorg$metadata
            cur.data <- reorg$exprs
        }

        # Counting across cells.
        out <- .recount_cells(exprs=cur.data, distance=distance, nsamples=nsamples, sample.id=sample.id-1L, 
                              cluster.centers=metadata$cluster.centers, cluster.info=metadata$cluster.info,
                              curcells=which(cell.id%%downsampling[r]==0L)-1L, markers=use.markers, filter=filter)
        if (is.character(out)) stop(out)
        out$cells <- .renew_indices(idx, out$cells)
        out$index <- idx[out$index]
        combined[[r]] <- out
    }        
    
    # Collating.
    combined.counts <- lapply(combined, "[[", name="counts")
    combined.coords <- lapply(combined, "[[", name="coords")
    combined.ass <- lapply(combined, "[[", name="cells")
    combined.center <- lapply(combined, "[[", name="index")
    combined.groups <- rep(seq_along(combined), lengths(combined.ass))

    combined.counts <- do.call(rbind, combined.counts)
    if (is.null(combined.counts)) combined.counts <- assay(x)[0,,drop=FALSE]
    combined.coords <- do.call(rbind, combined.coords)
    if (is.null(combined.coords)) combined.coords <- intensities(x)[0,used,drop=FALSE]
    combined.center <- unlist(combined.center)
    if (is.null(combined.center)) combined.center <- integer(0)
    combined.groups <- unlist(combined.groups)
    if (is.null(combined.groups)) combined.groups <- integer(0)

    # Returning a cyData object.
    cyData(assays=list(counts=combined.counts),
           intensities=combined.coords,
           cellAssignments=unlist(combined.ass, recursive=FALSE),
           cellIntensities=my.markers,
           markerData=markerData(x)[used,,drop=FALSE],
           cellData=cellData(x),
           rowData=DataFrame(center.cell=combined.center, group=combined.groups),
           colData=colData(x)
           ) 
}

.renew_indices <- function(ref, assign) {
    out <- .Call(cxx_renew_indices, ref, assign)
    if (is.character(out)) stop(out)
    return(out)
}
