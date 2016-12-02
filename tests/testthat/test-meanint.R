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
    to.use <- rep(c(TRUE, FALSE), each=5)
    filter <- 1L
    tol <- 0.5
    if (setup==2L) {
        to.use <- rep(c(TRUE, FALSE), c(8, 2))
    } else if (setup==3L) {
        to.use <- rep(c(TRUE, FALSE), each=5)
    } else if (setup==4L) {
        tol <- 1L
    }
    
    suppressWarnings(cd <- prepareCellData(fs, markers=to.use))
    out <- countCells(cd, filter=1L)
    rcnt <- meanIntensities(out, markers=!to.use)

    expect_identical(assay(out, "counts"), assay(rcnt, "counts"))
    expect_identical(cellData(out), cellData(rcnt))
    expect_identical(colData(out), colData(rcnt))
    expect_identical(rowData(out), rowData(rcnt))
    expect_identical(cellIntensities(out), cellIntensities(rcnt))
    expect_identical(markerData(rcnt), markerData(out))
    
    ref.groups <- unpackIndices(cellAssignments(out))
    for (u in which(!to.use)) { 
        cur.assay <- assay(rcnt, paste0("mean.", markernames(rcnt)[u]))
        ref.assay <- matrix(NA_real_, length(ref.groups), ncol(rcnt))
        colnames(ref.assay) <- sampleNames(rcnt)
        sample.names <- sampleNames(rcnt)[cellData(rcnt)$sample.id]

        for (r in seq_len(nrow(out))) {
            cur.group <- ref.groups[[r]]
            by.sample <- split(cur.group, sample.names[cur.group])

            for (s in names(by.sample)) { 
                ref.assay[r,s] <- mean(cellIntensities(rcnt)[u,by.sample[[s]]])
            }
        }
        expect_equal(ref.assay, cur.assay)
    } 
}



