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
        filter <- 5L
    } else if (setup==4L) {
        tol <- 1L
    }
    
    suppressWarnings(cd <- prepareCellData(fs, markers=to.use))
    out <- countCells(cd, filter=1L)
    ref.groups <- unpackIndices(cellAssignments(out))
    rcnt <- recountCells(out, markers=!to.use, filter=filter, tol=tol)

    expect_identical(cellData(out), cellData(rcnt))
    expect_identical(colData(out), colData(rcnt))
    expect_identical(cellIntensities(out), cellIntensities(rcnt))
    expect_identical(markerData(out), markerData(rcnt)[,1,drop=FALSE])
    expect_identical(!to.use, markerData(rcnt)[,2])
    expect_true(all(rowSums(assay(rcnt))>=filter))

    # Getting the cells currently in each nested hypersphere.
    all.ass <- unpackIndices(cellAssignments(out))
    collected <- list()
    for (r in seq_len(nrow(out))) {
        cur.group <- all.ass[[r]]
        central <- rowData(out)$center.cell[r]
 
        centre.point <- cellIntensities(out)[!to.use,central]
        cur.data <- cellIntensities(out)[!to.use,cur.group,drop=FALSE]
        keep <- sqrt(colSums((centre.point - cur.data)^2)) <= sqrt(sum(!to.use)) * tol
        collected[[r]] <- cur.group[keep]
    } 

    # Cross-referencing across hyperspheres.
    final.collected <- final.center <- final.group <- list()
    i <- 1L
    for (r in seq_len(nrow(out))) { 
        idx <- all.ass[[r]]
        all.centers <- intersect(idx, rowData(out)$center.cell)
        for (x in match(all.centers, rowData(out)$center.cell)) { 
            cur.collected <- intersect(idx, collected[[x]])
            if (length(cur.collected)< filter) next 
            final.collected[[i]] <- cur.collected
            final.center[[i]] <- rowData(out)$center.cell[x]
            final.group[[i]] <- r
            i <- i+1L
        }
    }

    # Matching up to the observed results.
    o <- match(paste0(unlist(final.group), ".", unlist(final.center)), 
               paste0(rowData(rcnt)$group, ".", rowData(rcnt)$center.cell))

    orcnt <- rcnt[o,]
    expect_identical(unpackIndices(cellAssignments(orcnt)), final.collected)
    expect_identical(rowData(orcnt)$center.cell, unlist(final.center))
    expect_identical(rowData(orcnt)$group, unlist(final.group))

    for (r in seq_len(nrow(orcnt))) {
        expect_identical(unname(assay(orcnt)[r,]), 
                         tabulate(cellData(orcnt)$sample.id[final.collected[[r]]], ncol(orcnt)))
    }
}



