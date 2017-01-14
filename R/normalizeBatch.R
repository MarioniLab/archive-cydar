normalizeBatch <- function(batch.x, batch.comp, mode=c("range", "quantile"), p=0.01)
# Performs quantile or range adjustment of different batches, given a 
# list of 'x' objects like that used for 'prepareCellData'
# and another list specifying the composition of samples per batch.
#
# written by Aaron Lun
# created 27 October 2016
# last modified 8 November 2016 
{
    if (is.null(batch.comp)) {
        batch.comp <- lapply(batch.x, function(i) rep(1, length(i)))
    }
    all.levels <- unique(unlist(batch.comp))
    batch.comp <- lapply(batch.comp, factor, levels=all.levels)
    nbatches <- length(batch.x)
    if (nbatches!=length(batch.comp)) {
        stop("length of 'batch.x' and 'batch.comp' must be identical")
    }

    # Computing the average number of samples from each batch to use in correction.
    comp.batches <- do.call(rbind, lapply(batch.comp, table))
    ref.comp <- colMeans(comp.batches)
    batch.weight <- ref.comp/comp.batches
    empty.factors <- colSums(!is.finite(batch.weight)) > 0
    if (all(empty.factors)) {
        stop("no level of 'batch.comp' is common to all batches")
    }
    use.batches <- colnames(batch.weight)[!empty.factors]

    # Checking the number of markers we're dealing with.
    batch.out <- list()
    for (b in seq_len(nbatches)) { 
        out <- .pull_out_data(batch.x[[b]])
        batch.out[[b]] <- out
        if (b==1L) {
            ref.markers <- out$markers
        } else if (!identical(ref.markers, out$markers)) { 
            stop("markers are not identical between batches")
        }
        if (length(out$samples)!=length(batch.comp[[b]])) {
            stop("corresponding elements of 'batch.comp' and 'batch.x' must have same lengths")
        }
    }

    # Computes sample- and batch-specific case weights to be used for all markers.
    batch.weights <- list()
    for (b in seq_len(nbatches)) { 
        cur.comp <- batch.comp[[b]]
        cur.out <- batch.out[[b]]
        cur.weights <- num.cells <- numeric(length(cur.out$exprs))
        
        for (s in seq_along(cur.out$exprs)) {
            sample.level <- as.character(cur.comp[s])
            num.cells[s] <- nrow(cur.out$exprs[[s]])
            if (sample.level %in% use.batches) { 
                cur.weights[s] <- 1/num.cells[s] * batch.weight[b,sample.level]
            }
        }

        batch.weights[[b]] <- rep(cur.weights, num.cells)
    }

    # Setting up an output object.
    output <- list()
    for (b in seq_len(nbatches)) { 
        cur.exprs <- list()   
        cur.out <- batch.out[[b]]
        for (s in seq_along(cur.out$samples)) {
            cur.exprs[[s]] <- cur.out$exprs[[s]]
            colnames(cur.exprs[[s]]) <- cur.out$markers
        }
        names(cur.exprs) <- cur.out$samples
        output[[b]] <- cur.exprs 
    }
    names(output) <- names(batch.x)

    mode <- match.arg(mode)
    for (m in ref.markers) {
        # Making a empirical cumprob -> quantile for each batch, given the marker.
        all.quanfun <- list()
        for (b in seq_len(nbatches)) { 
            cur.out <- batch.out[[b]]
            all.obs <- list()
            for (s in seq_along(cur.out$exprs)) {
                all.obs[[s]] <- cur.out$exprs[[s]][,m]
            }

            all.obs <- unlist(all.obs)
            o <- order(all.obs)
            all.obs <- all.obs[o]
            cur.weights <- batch.weights[[b]][o]

            # Taking the midpoint of each step, rather than the start/end points.
            mid.cum.weight <- cumsum(cur.weights) - cur.weights/2
            total.weight <- sum(cur.weights)
            all.quanfun[[b]] <- approxfun(mid.cum.weight/total.weight, all.obs, rule=2)
        }

        if (mode=="quantile") { 
            for (b in seq_len(nbatches)) { 
                current.quants <- environment(all.quanfun[[b]])$y
                current.probs <- environment(all.quanfun[[b]])$x

                # Correcting intensities for each sample in each batch; first by
                # computing the average distribution to which all others should be squeezed.
                cur.average <- 0
                for (b2 in seq_len(nbatches)) { 
                    if (b==b2) { 
                        cur.average <- cur.average + current.quants
                    } else{
                        cur.average <- cur.average + all.quanfun[[b2]](current.probs)
                    }
                }
                cur.average <- cur.average/length(batch.x)

                # Constructing a function to do that squeezing (we could use the ordering to 
                # get exactly which entry corresponds to which original observation, but that
                # doesn't allow for samples that weren't used, e.g., due to non-common levels).
                converter <- approxfun(current.quants, cur.average, rule=2) 
            
                cur.out <- batch.out[[b]]
                for (s in seq_along(cur.out$exprs)) {
                    output[[b]][[s]][,m] <- converter(cur.out$exprs[[s]][,m])                
                }
            }        
        } else {
            # Computing the average max/min.
            all.min <- all.max <- 0
            for (b in seq_len(nbatches)) { 
                all.min <- all.min + all.quanfun[[b]](p)
                all.max <- all.max + all.quanfun[[b]](1-p)
            }
            all.min <- all.min/nbatches
            all.max <- all.max/nbatches

            # Scaling intensities per batch so that the observed range equals the average range.
            for (b in seq_len(nbatches)) {
                targets <- c(all.min, all.max)
                current <- c(all.quanfun[[b]](p), all.quanfun[[b]](1-p))
                fit <- lm(targets ~ current)
                cur.out <- batch.out[[b]]

                for (s in seq_along(cur.out$exprs)) {
                    output[[b]][[s]][,m] <- cur.out$exprs[[s]][,m] * coef(fit)[2] + coef(fit)[1]
                }
            }
            
        }
    }
    return(output)
}

