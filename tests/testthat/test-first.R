# This checks that the findFirstSphere function is operating properly. 

require(cydar); require(testthat)
set.seed(400)
nhypers <- 1000
nmarkers <- 10
coords <- matrix(rnorm(nhypers*nmarkers, sd=1), nrow=nhypers, ncol=nmarkers)

for (threshold in c(0.5, 1, 2)) {
    keep <- findFirstSphere(coords, threshold=threshold)
    nkeep <- findFirstSphere(coords, threshold=threshold, naive=TRUE)
    expect_identical(keep, nkeep)

    tcoords <- t(coords)
    ref <- !logical(nhypers)
    keep <- logical(nhypers)
    for (i in seq_len(nhypers)) {
        if (i>1) {
            if (any(colSums(abs(tcoords[,keep,drop=FALSE] - tcoords[,i]) <= threshold)==nmarkers)) {
                ref[i] <- FALSE
            } else {
                keep[i] <- TRUE
            }
        } else {
            keep[i] <- TRUE
        }
    }

    expect_identical(ref, keep)
}

