normalizeBatch <- function(batch.x, batch.comp, mode=c("range", "warp"), p=0.01)
# Performs warp- or range-based adjustment of different batches, given a 
# list of 'x' objects like that used for 'prepareCellData'
# and another list specifying the composition of samples per batch.
#
# written by Aaron Lun
# created 27 October 2016
# last modified 30 January 2017
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
        # Putting together observations.
        all.obs <- list()
        for (b in seq_len(nbatches)) { 
            cur.out <- batch.out[[b]]
            cur.obs <- list()
            for (s in seq_along(cur.out$exprs)) {
                cur.obs[[s]] <- cur.out$exprs[[s]][,m]
            }
            all.obs[[b]] <- unlist(cur.obs)
        }

        if (mode=="warp") { 
            # Performing warp-based normalization to a reference.
            converters <- .transformDistr(all.obs, batch.weights, m)
            for (b in seq_len(nbatches)) { 
                converter <- converters[[b]]
                cur.out <- batch.out[[b]]
                for (s in seq_along(cur.out$exprs)) {
                    output[[b]][[s]][,m] <- converter(cur.out$exprs[[s]][,m])                
                }
            }        
        } else {
            # Computing the average max/min.
            batch.min <- batch.max <- numeric(nbatches)
            for (b in seq_len(nbatches)) { 
                cur.obs <- all.obs[[b]]
                o <- order(cur.obs)
                cur.obs <- cur.obs[o]
                cur.wts <- batch.weights[[b]][o]

                # Taking the midpoint of each step, rather than the start/end points.
                mid.cum.weight <- cumsum(cur.wts) - cur.wts/2
                total.weight <- sum(cur.wts)
                
                # Getting the left/right extreme.
                out <- approx(mid.cum.weight/total.weight, cur.obs, xout=c(p, 1-p), rule=2)$y
                batch.min[b] <- out[1]
                batch.max[b] <- out[2]
            }
            targets <- c(mean(batch.min), mean(batch.max))

            # Scaling intensities per batch so that the observed range equals the average range.
            for (b in seq_len(nbatches)) {
                current <- c(batch.min[b], batch.max[b])
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

.transformDistr <- function(all.obs, all.wts, name) {
    # Subsample intensities proportional to weights.
    cur.ffs <- list()
    for (b in seq_along(all.obs)) {
        cur.obs <- all.obs[[b]]
        cur.wts <- all.wts[[b]]
        chosen <- sample(cur.obs, length(cur.obs), prob=cur.wts, replace=TRUE)
        chosen <- c(chosen, range(cur.obs)) # Adding also the first and last entries.
        cur.ffs[[b]] <- flowFrame(cbind(M=chosen))
    }        

    names(cur.ffs) <- names(all.obs)
    fs <- as(cur.ffs, "flowSet")
    colnames(fs) <- name

    # Applying warping normalization, as described in the flowStats vignette.
    require(flowStats)
    norm <- normalization(normFunction=function(x, parameters, ...) { flowStats::warpSet(x, parameters, ...) },
                          parameters=name, arguments=list(monwrd=TRUE))
    new.fs <- normalize(fs, norm)

    # Defining warp functions (setting warpFuns doesn't really work, for some reason).
    converter <- list()
    for (b in seq_along(all.obs)) {
        old.i <- exprs(fs[[b]])[,1]
        new.i <- exprs(new.fs[[b]])[,1]
        converter[[b]] <- splinefun(old.i, new.i)
    }
    return(converter)
}

