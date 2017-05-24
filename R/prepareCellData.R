prepareCellData <- function(x, naive=FALSE, markers=NULL, ...) 
# Converts it into a format more suitable for high-speed analysis.
# Also does k-means clustering to generate the necessary clusters.
#
# written by Aaron Lun
# created 14 August 2016
# last modified 24 May 2017
{
    on.exit({gc()}) # Getting rid of huge memory structures that have built up.

    cell.data <- .pull_out_data(x)
    sample.names <-  cell.data$samples
    marker.names <- cell.data$markers
    exprs.list <- cell.data$exprs

    exprs <- do.call(rbind, exprs.list)
    sample.id <- rep(seq_along(exprs.list), sapply(exprs.list, nrow))
    cell.id <- unlist(lapply(exprs.list, function(x) { seq_len(nrow(x)) } ), use.names=FALSE)

    # Picking markers to use.
    used <- .chosen_markers(markers, marker.names)
    if (naive) { 
        exprs <- t(exprs)
        metadata <- list()
    } else {
        reorg <- .reorganize_cells(exprs, sample.id, cell.id, used, ...) 
        exprs <- reorg$exprs
        sample.id <- reorg$sample.id
        cell.id <- reorg$cell.id
        metadata <- reorg$metadata
    }
  
    # Collating the output. 
    CyData(cellIntensities=exprs,
           markerData=DataFrame(row.names=marker.names, used=used),
           colData=DataFrame(row.names=sample.names),
           assays=matrix(0L, 0, length(sample.names)),
           cellData=DataFrame(sample.id=sample.id, cell.id=cell.id),
           metadata=metadata)
}

.pull_out_data <- function(x)
# Pulling out data so we don't have to rely on ncdfFlowSet input.
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

.reorganize_cells <- function(exprs, sample.id, cell.id, used=NULL, ...) 
# Reorganizing for fast lookup via K-means clustering.
{
    if (is.null(used)) { 
        used <- rep(TRUE, ncol(exprs))
    }
    if (!all(used)) { 
        used.exprs <- exprs[,used,drop=FALSE]
    } else {
        used.exprs <- exprs
    }

    # Running K-means, with protection against the case where N=nrow.
    N <- ceiling(sqrt(nrow(exprs)))
    if (N==nrow(exprs)) { 
        out <- list(cluster=seq_len(N), centers=used.exprs)
    } else { 
        tryCatch({
            out <- suppressWarnings(kmeans(used.exprs, centers=N, ...))
        }, error=function(e) {
        }, finally={
            # Adding jitter, if many cells are duplicated.
            out <- suppressWarnings(kmeans(jitter(used.exprs), centers=N, ...))
        })
    }
    
    by.clust <- split(seq_len(nrow(exprs)), out$cluster)
    accumulated <- 0L
    nclust <- length(by.clust) # should be N, but redefining just in case...
    clust.info <- new.exprs <- new.samples <- new.cells <- vector("list", nclust)

    # Compiling to something that can be quickly accessed at the C++ level.
    for (clust in seq_len(nclust)) {
        chosen <- by.clust[[clust]]
        current.vals <- t(exprs[chosen,,drop=FALSE])
        cur.dist <- sqrt(colSums((out$centers[clust,] - current.vals[used,,drop=FALSE])^2))

        o <- order(cur.dist)
        new.exprs[[clust]] <- current.vals[,o,drop=FALSE]
        new.samples[[clust]] <- sample.id[chosen][o]
        new.cells[[clust]] <- cell.id[chosen][o]

        cur.dist <- cur.dist[o]
        clust.info[[clust]] <- list(accumulated, cur.dist)
        accumulated <- accumulated + length(o)
    }
   
    return(list(exprs=do.call(cbind, new.exprs), 
                metadata=list(cluster.centers=t(out$centers), cluster.info=clust.info),
                sample.id=unlist(new.samples),
                cell.id=unlist(new.cells)))
} 


