prepareCellData <- function(x, naive=FALSE, ...) 
# Converts it into a format more suitable for high-speed analysis.
# Also does k-means clustering to generate the necessary clusters.
#
# created by Aaron Lun
# created 14 August 2016
# last modified 14 November 2016    
{
    on.exit({gc()}) # Getting rid of huge memory structures that have built up.

    cell.data <- .pull_out_data(x)
    nsamples <- length(cell.data$samples)
    nmarkers <- length(cell.data$markers)
    exprs.list <- cell.data$exprs

    exprs <- do.call(rbind, exprs.list)
    sample.id <- rep(seq_along(exprs.list), sapply(exprs.list, nrow))
    cell.id <- unlist(lapply(exprs.list, function(x) { seq_len(nrow(x)) } ), use.names=FALSE)

if (!naive) {
    # K-means clustering.
    N <- ceiling(sqrt(nrow(exprs)))
    tryCatch({
        out <- kmeans(exprs, centers=N, ...)
    }, error=function(e) {
    }, finally={
        # Adding jitter, if many cells are duplicated.
        out <- kmeans(jitter(exprs), centers=N, ...)
    })
    by.clust <- split(seq_len(nrow(exprs)), out$cluster)
    clust.info <- list()
    new.exprs <- new.samples <- new.cells <- list()
    accumulated <- 0L

    # Compiling to something that can be quickly accessed at the C++ level.
    for (clust in seq_len(nrow(out$centers))) {
        chosen <- by.clust[[clust]]
        current.vals <- t(exprs[chosen,,drop=FALSE])
        cur.dist <- sqrt(colSums((out$centers[clust,] - current.vals)^2))

        o <- order(cur.dist)
        new.exprs[[clust]] <- current.vals[,o,drop=FALSE]
        new.samples[[clust]] <- sample.id[chosen][o]
        new.cells[[clust]] <- cell.id[chosen][o]

        cur.dist <- cur.dist[o]
        clust.info[[clust]] <- list(accumulated, cur.dist)
        accumulated <- accumulated + length(o)
    }

    all.exprs <- do.call(cbind, new.exprs)
    metadata <- list(cluster.centers=t(out$centers), cluster.info=clust.info)
    sample.id <- unlist(new.samples)
    cell.id <- unlist(new.cells)
} else {
    all.exprs <- t(exprs)
    metadata <- list()
} 
  
    # Collating the output. 
    cyData(cellIntensities=all.exprs,
           markerData=DataFrame(row.names=cell.data$markers),
           colData=DataFrame(row.names=cell.data$samples),
           assays=matrix(0L, 0, length(cell.data$samples)),
           cellData=DataFrame(sample.id=sample.id, cell.id=cell.id),
           metadata=metadata)
}

.check_cell_data <- function(x, check.clusters=TRUE) {
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
        stopifnot(nrow(cluster.centers)==nmarkers(x))
        
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

.pull_out_data <- function(x)
# Pulling out data so we don't have to rely on ncdfFlowSet.
{
    if (is.list(x)) { 
        sample.names <- names(x)
        if (is.null(sample.names)) { stop("list must be named by sample") }
        marker.names <- colnames(x[[1]])
        for (i in sample.names) {
            x[[i]] <- as.matrix(x[[i]])
            tmp <- colnames(x[[i]])
            if (is.null(tmp)) { stop("column names must be labelled with marker identities"); }
            stopifnot(identical(tmp, marker.names))
        }
        expr.val <- x
    } else if (is(x, "ncdfFlowSet")) { 
        sample.names <- Biobase::sampleNames(x)
        marker.names <- BiocGenerics::colnames(x)
        by.sample <- seq_along(sample.names)
        expr.val <- lapply(by.sample, FUN=function(i) flowCore::exprs(x[[i]]))
    } else {
        stop("'cell.data' must be a list or ncdfFlowSet object") 
    }
    return(list(samples=sample.names, markers=marker.names, exprs=expr.val))
}

