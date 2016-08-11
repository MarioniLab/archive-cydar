countCells <- function(cell.data, tol=0.5, BPPARAM=bpparam(), downsample=10, filter=10, naive=FALSE)
# Extracts counts for each cell in a CyTOF experiment, based on the number of surrounding cells 
# from each sample, given a prepared set of expression values for all cells in each sample.
#
# written by Aaron Lun
# created 21 April 2016
# last modified 11 August 2016
{
    .check_cell_data(cell.data)
    sample.id <- attributes(cell.data)$sample.id
    cell.id <- attributes(cell.data)$cell.id
    if (naive) {
        cluster.centers <- cluster.info <- NULL
    } else {
        cluster.centers <- attributes(cell.data)$cluster.centers
        cluster.info <- attributes(cell.data)$cluster.info
    }
    samples <- attributes(cell.data)$samples

    # Scaling the distance by the number of parameters. 
    markers <- attributes(cell.data)$markers
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
    out <- bplapply(allocations, FUN=.recount_cells, exprs=cell.data, nsamples=length(samples), 
                    sample.id=sample.id, distance=distance, cluster.centers=cluster.centers, 
                    cluster.info=cluster.info, filter=filter, BPPARAM=BPPARAM)

    # Assembling output into a coherent set.
    out.counts <- out.coords <- list()
    for (i in seq_along(out)) {
        if (is.character(out[[i]])) { stop(out[[i]]) }
        out.counts[[i]] <- out[[i]]$counts
        out.coords[[i]] <- out[[i]]$coords
    }

    out.counts <- do.call(rbind, out.counts)
    colnames(out.counts) <- samples
    out.coords <- do.call(rbind, out.coords)
    colnames(out.coords) <- markers

    # Ordering them (not strictly necessary, just for historical reasons).
    i <- as.integer(sub("^c", "", rownames(out.counts))) 
    o <- order(sample.id[i], cell.id[i])
    out.counts <- out.counts[o,,drop=FALSE]
    out.coords <- out.coords[o,,drop=FALSE]

    all.ncells <- tabulate(sample.id+1L, length(samples))
    return(list(counts=out.counts, coordinates=out.coords, totals=all.ncells, tolerance=tol))
}

.recount_cells <- function(curcells, exprs, nsamples, sample.id, distance, cluster.centers, cluster.info, filter) 
# Helper function so that BiocParallel call is self-contained.
{
    out <- .Call(cxx_count_cells, exprs, distance, nsamples, sample.id, cluster.centers, cluster.info, curcells)
    if (!is.character(out)) { 
        counts <- out[[1]]
        coords <- out[[2]]
        keep <- rowSums(counts) >= filter
        counts <- counts[keep,,drop=FALSE]
        coords <- coords[keep,,drop=FALSE]
        if (any(keep)) { rownames(counts) <- rownames(coords) <- paste0("c", curcells[keep] + 1L) }
        return(list(counts=counts, coords=coords))
    }
    return(out)
}

