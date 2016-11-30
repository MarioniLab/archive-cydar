#####################################################
# This tests the count machinery, to make sure the counts are right.

require(cydar); require(testthat)

set.seed(100)
nmarkers <- 10
tol <- 0.5

# Setup.
ncells1 <- 1001
all.values1 <- matrix(rnorm(ncells1*nmarkers, sd=1), nrow=ncells1, ncol=nmarkers)
colnames(all.values1) <- paste0("X", seq_len(nmarkers))

ncells2 <- 2001
all.values2 <- matrix(rnorm(ncells2*nmarkers, sd=1), nrow=ncells2, ncol=nmarkers)
colnames(all.values2) <- colnames(all.values1)
fs <- list(A=all.values1, B=all.values2)

# Counting with specified parameters.

for (setup in 1:5) { 
    to.use <- rep(c(TRUE, FALSE), each=5)
    filter <- 1L
    downsample <- NULL
    tol <- 0.5
    if (setup==2L) {
        to.use <- rep(c(TRUE, FALSE), c(8, 2))
    } else if (setup==3L) {
        filter <- 5L
    } else if (setup==4L) {
        downsample <- 10L
    } else if (setup==5L) {
        tol <- 1L
    }
    
    suppressWarnings(cd <- prepareCellData(fs, markers=to.use))
    out <- countCells(cd, filter=1L)
    ref.groups <- unpackIndices(cellAssignments(out))
    rcnt <- recountCells(out, markers=!to.use, filter=filter, downsample=downsample, tol=tol)

    expect_identical(cellData(out), cellData(rcnt))
    expect_identical(colData(out), colData(rcnt))
    expect_identical(cellIntensities(out), cellIntensities(rcnt))
    expect_identical(markerData(out), markerData(rcnt)[,1,drop=FALSE])
    expect_identical(!to.use, markerData(rcnt)[,2])
    expect_true(all(rowSums(assay(rcnt))>=filter))

    for (r in seq_len(nrow(rcnt))) {
        idx <- unpackIndices(cellAssignments(rcnt)[r])[[1]]
        central <- rowData(rcnt)$center.cell[r]

        # Checking cells are in the specified group.
        cur.group <- ref.groups[[rowData(rcnt)$group[r]]]
        expect_true(central %in% cur.group)
        expect_true(central %in% idx)
        expect_true(all(idx %in% cur.group))

        # Checking currently assigned cells are correct.
        centre.point <- cellIntensities(out)[!to.use,central]
        cur.data <- cellIntensities(out)[!to.use,cur.group,drop=FALSE]
        keep <- sqrt(colSums((centre.point - cur.data)^2)) <= sqrt(sum(!to.use)) * tol
        expect_identical(cur.group[keep], idx)

        # Checking that they yield the specified counts.
        expect_equivalent(assay(rcnt)[r,], tabulate(cellData(out)$sample.id[idx], nbin=ncol(out)))
    }
}



