prepareCellData <- function(x, ...) 
# Converts it into a format more suitable for high-speed analysis.
# Also does k-means clustering to generate the necessary clusters.
{
    on.exit({gc()}) # Getting rid of huge memory structures that have built up.

    cell.data <- .pull_out_data(x)
    nsamples <- length(cell.data$samples)
    nmarkers <- length(cell.data$markers)
    exprs.list <- cell.data$exprs

    exprs <- do.call(rbind, exprs.list)
    sample.id <- rep(seq_along(exprs.list) - 1L, sapply(exprs.list, nrow))
    cell.id <- unlist(lapply(exprs.list, function(x) { seq_len(nrow(x)) } )) - 1L

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
  
    # Collating the output. 
    output <- do.call(cbind, new.exprs)
    attributes(output)$markers <- cell.data$markers
    attributes(output)$samples <- cell.data$samples
    attributes(output)$sample.id <- unlist(new.samples)
    attributes(output)$cell.id <- unlist(new.cells)
    attributes(output)$cluster.centers <- t(out$centers)
    attributes(output)$cluster.info <- clust.info
    return(output)
}

.check_cell_data <- function(cell.data) {
    markers <- attributes(cell.data)$markers
    stopifnot(length(markers)==nrow(cell.data))

    samples <- attributes(cell.data)$samples
    stopifnot(!is.null(samples))
    
    sample.id <- attributes(cell.data)$sample.id
    stopifnot(all(sample.id >= 0L & sample.id < length(samples)))
    stopifnot(length(sample.id)==ncol(cell.data))

    cell.id <- attributes(cell.data)$cell.id
    stopifnot(length(cell.id)==ncol(cell.data))

    cluster.centers <- attributes(cell.data)$cluster.centers        
    stopifnot(nrow(cluster.centers)==nrow(cell.data))
    
    cluster.info <- attributes(cell.data)$cluster.info
    stopifnot(ncol(cluster.centers)==length(cluster.info))
    for (x in cluster.info) {
        stopifnot(!is.unsorted(x[[2]]))
        stopifnot(x[[1]] >= 0L & x[[1]]+length(x[[2]]) <= ncol(cell.data))
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

