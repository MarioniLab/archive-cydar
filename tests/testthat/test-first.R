# This checks that the findFirstSphere function is operating properly. 

require(cydar); require(testthat)
set.seed(400)
nhypers <- 1000

for (nmarkers in c(10, 20, 30)) {
    coords <- matrix(rnorm(nhypers*nmarkers, sd=1), nrow=nhypers, ncol=nmarkers)
    pval <- rbeta(nhypers, 1, 10)
    tcoords <- t(coords)
    o <- order(pval) 
    block <- sample(3, nhypers, replace=TRUE)
    
    for (threshold in c(0.5, 1, 2)) {
        xkeep <- findFirstSphere(coords, pvalues=pval, threshold=threshold)
        nkeep <- findFirstSphere(coords, pvalues=pval, threshold=threshold, naive=TRUE)
        expect_identical(xkeep, nkeep)
        
        ref <- logical(nhypers)
        locals <- integer(nhypers)
        for (j in seq_along(o)) {
            i <- o[j]
            if (j > 1) {
                if (!any(colSums(abs(tcoords[,ref,drop=FALSE] - tcoords[,i]) <= threshold)==nmarkers)) {
                    ref[i] <- TRUE
                }
            } else {
                ref[i] <- TRUE
            }
        }
        
        expect_identical(ref, xkeep)

        # Trying with a blocking factor.
        bkeep <- findFirstSphere(coords, pvalues=pval, threshold=threshold, block=block)
        all.out <- logical(nhypers)
        for (i in unique(block)) {
            chosen <- block==i
            all.out[chosen] <- findFirstSphere(coords[chosen,,drop=FALSE], pvalues=pval[chosen], threshold=threshold)
        }
        expect_identical(all.out, bkeep)
    }
}

