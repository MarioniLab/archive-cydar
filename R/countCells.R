countCells <- function(x, tol=0.5, BPPARAM=bpparam(), downsample=10, filter=10, naive=FALSE)
# Extracts counts for each cell in a CyTOF experiment, based on the number of surrounding cells 
# from each sample, given a prepared set of expression values for all cells in each sample.
#
# written by Aaron Lun
# created 21 April 2016
# last modified 14 November 2016
{
    .check_cell_data(x, check.clusters=!naive)
    sample.id <- cellData(x)$sample.id - 1L # Get to zero indexing.
    cell.id <- cellData(x)$cell.id - 1L
    if (naive) {
        cluster.centers <- cluster.info <- NULL
    } else {
        cluster.centers <- metadata(x)$cluster.centers
        cluster.info <- metadata(x)$cluster.info
    }
    samples <- colnames(x)

    # Scaling the distance by the number of parameters. 
    markers <- rownames(markerData(x))
    nmarkers <- length(markers)
    distance <- tol * sqrt(nmarkers) 
    if (distance <= 0) {
        warning("setting a non-positive distance to a small offset")
        distance <- 1e-8        
    }
    
    # Getting some other values.
    downsample <- as.integer(downsample)
    chosen <- which((cell.id %% downsample) == 0L) - 1L
    N <- ceiling(length(chosen)/bpworkers(BPPARAM))
    core.assign <- rep(seq_len(bpworkers(BPPARAM)), each=N, length.out=length(chosen))
    allocations <- split(chosen, core.assign)

    # Parallel analysis.
    ci <- cellIntensities(x)
    out <- bplapply(allocations, FUN=.recount_cells, exprs=ci, nsamples=length(samples), 
                    sample.id=sample.id, distance=distance, cluster.centers=cluster.centers, 
                    cluster.info=cluster.info, filter=filter, BPPARAM=BPPARAM)

    # Assembling output into a coherent set.
    out.counts <- out.coords <- out.cells <- out.index <- list()
    for (i in seq_along(out)) {
        if (is.character(out[[i]])) { stop(out[[i]]) }
        out.counts[[i]] <- out[[i]]$counts
        out.coords[[i]] <- out[[i]]$coords
        out.cells[[i]] <- out[[i]]$cells
        out.index[[i]] <- out[[i]]$index
    }

    out.counts <- do.call(rbind, out.counts)
    colnames(out.counts) <- samples
    out.coords <- do.call(rbind, out.coords)
    colnames(out.coords) <- markers
    out.cells <- unlist(out.cells, recursive=FALSE)
    names(out.cells) <- NULL
    out.index <- unlist(out.index)
    names(out.index) <- NULL

    # Ordering them (not strictly necessary, just for historical reasons).
    o <- order(sample.id[out.index], cell.id[out.index])
    out.counts <- out.counts[o,,drop=FALSE]
    rownames(out.counts) <- seq_along(o)
    out.coords <- out.coords[o,,drop=FALSE]
    out.cells <- out.cells[o]
    out.index <- out.index[o]

    all.ncells <- tabulate(sample.id+1L, length(samples))
    output <- new("cyData", x, assays=Assays(list(counts=out.counts)), 
                  intensities=out.coords, cellAssignments=out.cells,
                  elementMetadata=DataFrame(center.cell=out.index)) 
    output$totals <- all.ncells
    metadata(output)$tol <- tol
    return(output)
}

.recount_cells <- function(curcells, exprs, nsamples, sample.id, distance, cluster.centers, cluster.info, filter) 
# Helper function so that BiocParallel call is self-contained.
{
    out <- .Call(cxx_count_cells, exprs, distance, nsamples, sample.id, cluster.centers, cluster.info, curcells)
    if (!is.character(out)) { 
        counts <- out[[1]]
        coords <- out[[2]]
        cells <- out[[3]]
        keep <- rowSums(counts) >= filter
        counts <- counts[keep,,drop=FALSE]
        coords <- coords[keep,,drop=FALSE]
        cells <- cells[keep]
        index <- curcells[keep] + 1L
        return(list(counts=counts, coords=coords, cells=cells, index=index))
    }
    return(out)
}

