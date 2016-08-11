#####################################################
# This tests the count machinery, to make sure the counts are right.

require(cyder); require(testthat)

set.seed(100)
for (setup in 1:5) {
    # Running our function.
    nmarkers <- 10
    downsample <- 10L
    downsample <- 10L
    tol <- 0.5
    filter <- 1L
    if (setup==2L) {
        nmarkers <- 5L
    } else if (setup==3L) {
        downsample <- 3L
    } else if (setup==4L) {
        tol <- 1L
    } else if (setup==5L) {
        tol <- 1L
        filter <- 5L
    }

    # Setup.
    ncells1 <- 1000
    all.values1 <- matrix(rnorm(ncells1*nmarkers, sd=1), nrow=ncells1, ncol=nmarkers)
    colnames(all.values1) <- paste0("X", seq_len(nmarkers))
    
    ncells2 <- 2000
    all.values2 <- matrix(rnorm(ncells2*nmarkers, sd=1), nrow=ncells2, ncol=nmarkers)
    colnames(all.values2) <- colnames(all.values1)
    fs <- list(A=all.values1, B=all.values2)

    # Counting with specified parameters.
    suppressWarnings(cd <- prepareCellData(fs))
    out <- countCells(cd, filter=filter, downsample=downsample, tol=tol)
    out2 <- countCells(cd, filter=filter, downsample=downsample, tol=tol, naive=TRUE)
    expect_identical(out, out2)

    # Reference counter.
    to.select1 <- seq_len(nrow(all.values1)) %% downsample == 1L
    to.select2 <- seq_len(nrow(all.values2)) %% downsample == 1L
    to.select <- c(to.select1, to.select2)
    combined <- rbind(all.values1, all.values2)
    origin <- rep(1:2, c(nrow(all.values1), nrow(all.values2)))
    new.dist <- tol * sqrt(nmarkers)

    collected.counts <- list()
    collected.meds <- list()
    index <- 1L
    for (i in which(to.select)) { 
        curdist <- sqrt(colSums((t(combined) - combined[i,])^2))
        inrange <- curdist <= new.dist
        collected.counts[[index]] <- tabulate(origin[inrange], nbins=length(fs))
        collected.meds[[index]] <- apply(combined, 2, FUN=function(x) { median(x[inrange]) })
        index <- index + 1L
    }
    collected.counts <- do.call(rbind, collected.counts)
    keep <- rowSums(collected.counts) >= filter
    
    # Comparison.
    colnames(collected.counts) <- names(fs)
    rownames(collected.counts) <- match(c(paste0("0.", which(to.select1)-1L), paste0("1.", which(to.select2)-1L)),
                                        paste0(attributes(cd)$sample.id, ".", attributes(cd)$cell.id))
    expect_identical(collected.counts[keep,,drop=FALSE], out$counts)

    collected.meds <- do.call(rbind, collected.meds)
    colnames(collected.meds) <- colnames(all.values1)
    rownames(collected.meds) <- rownames(collected.counts)
    expect_equal(collected.meds[keep,,drop=FALSE], out$coordinates)

    expect_equal(tol, out$tolerance)
    expect_equal(c(ncells1, ncells2), out$totals)
}

# Running silly settings.
suppressWarnings(out <- countCells(cd, filter=1L, tol=0))
expect_true(all(rowSums(out$counts)==1L))
out <- countCells(cd, filter=Inf)
expect_identical(nrow(out$counts), 0L)
