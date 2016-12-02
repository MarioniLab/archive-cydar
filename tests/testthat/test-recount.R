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

for (setup in 1:4) { 
    to.use.1 <- rep(c(TRUE, FALSE), each=5)
    to.use.2 <- !to.use.1
    filter <- 1L
    tol <- 0.5
    if (setup==2L) {
        to.use.1 <- rep(c(TRUE, FALSE), c(8, 2))
        to.use.2 <- !to.use.1
    } else if (setup==3L) {
        to.use.1 <- rep(c(TRUE, FALSE), each=5)
        to.use.2 <- rep(c(FALSE, TRUE, FALSE), c(5, 2, 3))
    } else if (setup==4L) {
        tol <- 1L
    }
    
    suppressWarnings(cd <- prepareCellData(fs, markers=to.use.1))
    out <- countCells(cd, filter=1L)
    rcnt <- recountCells(out, markers=to.use.2, tol=tol)

    expect_identical(cellData(out), cellData(rcnt))
    expect_identical(colData(out), colData(rcnt))
    expect_identical(rowData(out), rowData(rcnt))
    expect_identical(cellIntensities(out), cellIntensities(rcnt))
    to.use <- to.use.1| to.use.2
    expect_identical(markerData(rcnt)$used, to.use)

    # Getting the cells currently in each nested hypersphere.
    all.ass <- unpackIndices(cellAssignments(out))
    collected <- list()
    for (r in seq_len(nrow(out))) {
        cur.group <- all.ass[[r]]
        central <- rowData(out)$center.cell[r]
 
        centre.point <- cellIntensities(out)[to.use,central]
        cur.data <- cellIntensities(out)[to.use,cur.group,drop=FALSE]
        keep <- sqrt(colSums((centre.point - cur.data)^2)) <= sqrt(sum(to.use)) * metadata(rcnt)$tol
        collected[[r]] <- cur.group[keep]
    } 

    expect_identical(unpackIndices(cellAssignments(rcnt)), collected)
    for (r in seq_len(nrow(rcnt))) {
        expect_identical(unname(assay(rcnt)[r,]), 
                         tabulate(cellData(rcnt)$sample.id[collected[[r]]], ncol(rcnt)))
    }
}



