diffIntDistr <- function(..., markers=NULL, npts=200) 
# Computes several matrices of differences between samples
# in the same and different batches. Distances in the intensity
# distributions for each marker are computed as KS statistics.
#
# written by Aaron Lun
# created 13 January 2017
# last modified 21 March 2017 
{
    # Setting up inputs.
    all.out <- list(...)
    if (is.null(names(all.out))) { names(all.out) <- "" }
    names.absent <- is.na(names(all.out)) | names(all.out)==""
    names(all.out)[names.absent] <-  paste0("Batch", which(names.absent))

    batch.out <- collected.markers <- sample.names <- vector("list", length(all.out))
    for (b in seq_along(all.out)) {
        current <- .pull_out_data(all.out[[b]])
        batch.out[[b]] <- current$exprs
        collected.markers[[b]] <- current$markers
        sample.names[[b]] <- paste0(names(all.out)[b], ".", current$samples)
    }
    
    batch.out <- unlist(batch.out, recursive=FALSE)
    sample.names <- unlist(sample.names)
    if (is.null(markers)) { 
        markers <- Reduce(intersect, collected.markers)
    }
   
    output <- vector("list", length(markers))
    names(output) <- markers
    for (m in markers) {  
        # Getting the min/max range for the current marker.
        curmin <- curmax <- NULL
        for (s in seq_along(batch.out)) {
            curmin <- min(curmin, batch.out[[s]][,m]) 
            curmax <- max(curmax, batch.out[[s]][,m]) 
        }

        # Computing densities for all samples.
        pts <- seq(from=curmin - 1e-6, to=curmax + 1e-6, length.out=npts)
        densities <- vector("list", length(batch.out))
        for (s in seq_along(batch.out)) {
            current <- batch.out[[s]][,m]
            densities[[s]] <- diff(findInterval(pts, sort(current)))/length(current)
        }

        # Computing sum of absolute differences between all pairs of samples.
        intermat <- matrix(0, length(batch.out), length(batch.out))
        for (s1 in seq_along(batch.out)) {
            for (s2 in seq_len(s1-1L)) {
                stat <- sum(abs(densities[[s1]] - densities[[s2]]))/2
                intermat[s1,s2] <- intermat[s2,s1] <- stat
            }
        }
        
        intermat <- intermat * diff(pts)[1]    
        rownames(intermat) <- colnames(intermat) <- sample.names
        output[[m]] <- intermat
    }

    # Setting up an indicator matrix to specify which goes where.
    ind <- matrix(0L, length(batch.out), length(batch.out))
    last.pos <- 0L
    for (b in seq_along(all.out)) {
        curx <- seq_along(all.out[[b]]) + last.pos
        ind[curx, curx] <- b
        last.pos <- last.pos + length(all.out[[b]])
    }
    rownames(ind) <- colnames(ind) <- sample.names

    return(list(index=ind, difference=output))
}

