#####################################################
# This tests the count machinery, to make sure the counts are right.

require(cydar); require(testthat)

set.seed(100)
for (setup in 1:5) {
    # Running our function.
    nmarkers <- 10
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
    ncells1 <- 1001
    all.values1 <- matrix(rnorm(ncells1*nmarkers, sd=1), nrow=ncells1, ncol=nmarkers)
    colnames(all.values1) <- paste0("X", seq_len(nmarkers))
    
    ncells2 <- 2001
    all.values2 <- matrix(rnorm(ncells2*nmarkers, sd=1), nrow=ncells2, ncol=nmarkers)
    colnames(all.values2) <- colnames(all.values1)
    fs <- list(A=all.values1, B=all.values2)

    # Counting with specified parameters.
    suppressWarnings(cd <- prepareCellData(fs))
    out <- countCells(cd, filter=filter, downsample=downsample, tol=tol)
    out2 <- countCells(cd, filter=filter, downsample=downsample, tol=tol, naive=TRUE)
    expect_equal(out, out2)

    # Reference counter.
    to.select1 <- seq_len(nrow(all.values1)) %% downsample == 1L
    to.select2 <- seq_len(nrow(all.values2)) %% downsample == 1L
    to.select <- c(to.select1, to.select2)
    combined <- rbind(all.values1, all.values2)
    origin <- rep(1:2, c(nrow(all.values1), nrow(all.values2)))
    cell.id <- c(seq_len(nrow(all.values1)), seq_len(nrow(all.values2)))
    new.dist <- tol * sqrt(nmarkers)

    collected.counts <- list()
    collected.meds.upper <- collected.meds.lower <- list()
    index <- 1L
    for (i in which(to.select)) { 
        curdist <- sqrt(colSums((t(combined) - combined[i,])^2))
        inrange <- curdist <= new.dist
        cursamples <- origin[inrange]
        collected.counts[[index]] <- tabulate(cursamples, nbins=length(fs))

        # Need to calculate lower/upper bounds for the median, as numerical imprecision has big effects.
        cur.meds <- apply(combined, 2, 
                          FUN=function(x) { 
                              x <- x[inrange]
                              w <- 1/c(ncells1, ncells2)[origin[inrange]]
                              o <- order(x)
                              x <- x[o]
                              w <- w[o]
                              p <- cumsum(w)/sum(w)
                              x[c(sum(p < 0.499999), sum(p < 0.500001))+1]
                          })
        collected.meds.lower[[index]] <- cur.meds[1,]
        collected.meds.upper[[index]] <- cur.meds[2,]

        index <- index + 1L
    }
    collected.counts <- do.call(rbind, collected.counts)
    keep <- rowSums(collected.counts) >= filter
    
    # Comparison.
    obs <- assay(out)
    dimnames(obs) <- NULL
    expect_identical(collected.counts[keep,,drop=FALSE], obs)
    expect_identical(rowData(out)$center.cell, match(c(paste0("1.", which(to.select1)), paste0("2.", which(to.select2)))[keep],
        paste0(cellData(cd)$sample.id, ".", cellData(cd)$cell.id)))

    med.coords <- intensities(out)
    collected.meds.lower <- do.call(rbind, collected.meds.lower)[keep,,drop=FALSE]
    collected.meds.upper <- do.call(rbind, collected.meds.upper)[keep,,drop=FALSE]
    expect_identical(dim(collected.meds.lower), dim(med.coords))
    expect_identical(dim(collected.meds.upper), dim(med.coords))
    out.of.range <- collected.meds.lower > med.coords | collected.meds.upper < med.coords
    expect_true(!any(out.of.range))

    expect_equal(tol, metadata(out)$tol)
    expect_equal(c(ncells1, ncells2), out$totals)

    # Checking the consistency of the compressed vectors with the counts.
    for (r in seq_along(cellAssignments(out))) { 
        cell.indices <- unpackIndices(cellAssignments(out)[r])
        expect_equivalent(assay(out)[r,], tabulate(cellData(out)$sample.id[cell.indices[[1]]], nbin=2))
        expect_identical(packIndices(cell.indices), cellAssignments(out)[r])
    }

    # Checking the consistency of the counts if only a subset of markers are used.
    chosen.markers <- c(2, 3, 4)
    fs.sub <- fs
    for (f in seq_along(fs.sub)) { fs.sub[[f]] <- fs.sub[[f]][,chosen.markers,drop=FALSE] }
    suppressWarnings(cd.sub <- prepareCellData(fs.sub))
    out.sub <- countCells(cd.sub, filter=filter, downsample=downsample, tol=tol)
    suppressWarnings(cd.sub.2 <- prepareCellData(fs, markers=chosen.markers))
    out.sub.2 <- countCells(cd.sub.2, filter=filter, downsample=downsample, tol=tol)
    expect_identical(assay(out.sub), assay(out.sub.2))
    for (r in seq_along(cellAssignments(out.sub))) {
        ix1 <- cellData(out.sub)[unpackIndices(cellAssignments(out.sub)[r])[[1]],]
        ix2 <- cellData(out.sub.2)[unpackIndices(cellAssignments(out.sub.2)[r])[[1]],]
        ix1 <- ix1[order(ix1$sample.id, ix1$cell.id),]
        ix2 <- ix2[order(ix2$sample.id, ix2$cell.id),]
        expect_identical(ix1, ix2)
    }
}

# Running silly settings.
suppressWarnings(out <- countCells(cd, filter=1L, tol=0))
expect_true(all(rowSums(assay(out))==1L))
out <- countCells(cd, filter=Inf)
expect_identical(nrow(out), 0L)

