normalizeBatch <- function(batch.x, batch.comp, mode=c("range", "curve"),
    p=0.01, ref=1, npts=512, extra.init=NULL)
# Performs curve- or range-based adjustment of different batches, given a 
# list of 'x' objects like that used for 'prepareCellData'
# and another list specifying the composition of samples per batch.
#
# written by Aaron Lun
# created 27 October 2016
# last modified 27 January 2017
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
        print(m)
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

        if (mode=="curve") { 
            # Performing curve-based normalization to a reference.
            converters <- .transformDistr(all.obs, batch.weights, ref=ref, npts=npts, extra.init=extra.init)
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
            all.min <- mean(batch.min)
            all.max <- mean(batch.max)

            # Scaling intensities per batch so that the observed range equals the average range.
            for (b in seq_len(nbatches)) {
                targets <- c(all.min, all.max)
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

.createCompareDistrFun <- function(ref.obs, cur.obs, ref.wts, cur.wts, npts=512) {
    x <- list(ref.obs, cur.obs)
    w <- list(ref.wts, cur.wts)
    npts <- as.integer(npts)

    # Handling zeroes separately.
    sum.zero <- sum.nzero <- numeric(length(x))
    for (b in seq_along(x)) { 
        cur.o <- x[[b]]
        cur.w <- w[[b]]
        is.zero <- cur.o < 1e-16
        sum.zero[b] <- sum(cur.w[is.zero])
        x[[b]] <- cur.o[!is.zero]
        cur.w <- cur.w[!is.zero]
        w[[b]] <- cur.w/sum(cur.w)
        sum.nzero[b] <- sum(w[[b]])
    }
    total.wts <- sum.zero + sum.nzero
    prop.zero <- sum.zero/total.wts
    prop.nzero <- sum.nzero/total.wts

    # Computing difference in PDFs (function optimized for speed w.r.t. with zeros).
    ref.obs <- x[[1]]
    ref.max <- max(ref.obs)
    ref.wts <- w[[1]]
    cur.obs <- x[[2]]
    cur.wts <- w[[2]]

    FUN <- function(new.obs) {
        maxpt <- max(ref.max, max(new.obs)) + 1
        pt.size <- maxpt/(npts-1)
        
        is.zero <- new.obs < 1e-16
        mod <- 0
        if (any(is.zero)) {
            mod <- sum(is.zero)/total.wts[2]
            new.obs <- new.obs[!is.zero]
            cur.wts <- cur.wts[!is.zero]
            cur.wts <- cur.wts/sum(cur.wts)
        }

        rout <- safe_density(x=ref.obs, weights=ref.wts, from=0, to=maxpt, n=npts)
        nout <- safe_density(x=new.obs, weights=cur.wts, from=0, to=maxpt, n=npts)
        sum(abs(rout$y * prop.nzero[1] - nout$y * (prop.nzero[2] - mod))) * pt.size + abs(prop.zero[1] - (prop.zero[2] + mod))
    }

    # Force intensities to be non-negative:
    return(list(FUN=FUN, rest.obs=pmax(cur.obs, 0))) 
}

safe_density <- function(x, ..., n) {
    # Protect against errors in edge cases.
    if (length(x)>=2L) {
        return(density(x=x, ...))
    } else if (length(x)==1L) {
        return(density(x=x, bw=0.1, ...))
    } else {
        return(list(y=numeric(n)))
    }
}

.transformDistr <- function(all.obs, all.wts, ref=1, npts=npts, extra.init=extra.init) {
    # Adjustment function is: y = A/(1+exp(B-exp(C)*x)) - A/(1+exp(B))
    # B and C are in the exponents so that exp(B) and exp(C) are always positive.
    transFUN <- function(par, x) {
        adj <- par[["A"]]/(1 + exp(par[["B"]]-exp(par[["C"]])*x)) - par[["A"]]/(1 + exp(par[["B"]]))
        pmax(0, x + adj)
    }

    # Figuring if the reference.
    ref <- as.integer(ref)
    if (ref <= 0L || ref >= length(all.obs)) {
        stop("'ref' must be one of the batches")
    }
    ref.obs <- all.obs[[ref]]
    ref.wts <- all.wts[[ref]]
    
    # Normalizing everything to the reference batch.
    converter <- list()
    for (b in seq_along(all.obs)) {
        if (b==ref) { 
            converter[[b]] <- function(x) x
            next 
        }
        cur.obs <- all.obs[[b]]
        cur.wts <- all.wts[[b]]

        # Generating the function to optimize.
        comp.out <- .createCompareDistrFun(ref.obs, cur.obs, ref.wts, cur.wts, npts=npts)
        compFUN <- comp.out$FUN
        rest.obs <- comp.out$rest.obs
                          
        # Running optim from many starting points (different shapes), to explore the lumpy space.
        last.val <- Inf
        opt.out <- NULL
        for (init in list(c(A=0,B=1,C=1), 
                          c(A=0,B=10,C=1),
                          c(A=0,B=1,C=-10),
                          extra.init)) { 
            cur.out <- optim(init, fn=function(par) {
                new.obs <- transFUN(par, rest.obs)
                compFUN(new.obs)
            })
            if (cur.out$value < last.val) {
                last.val <- cur.out$value
                opt.out <- cur.out
            }
        }

        # Slight contortion to avoid lexical scoping.
        print(opt.out$par)
        converter[[b]] <- (function(par) { par; function(x) transFUN(par, x) })(opt.out$par)
    }
    
    return(converter)
}
