#####################################################
# Testing the plotting functions.

require(cydar); require(testthat)
set.seed(200)

x <- runif(1000)
y <- runif(1000)
lfc <- rnorm(1000)

out <- plotCellLogFC(x, y, lfc)
maxchange <- max(abs(lfc))
expect_equal(maxchange, -as.double(head(names(out), 1)))
expect_equal(maxchange, as.double(tail(names(out), 1)))

out <- plotCellLogFC(x, y, lfc, max.logFC=2)
maxchange <- 2
expect_equal(maxchange, -as.double(head(names(out), 1)))
expect_equal(maxchange, as.double(tail(names(out), 1)))

# Same for plotting on intensity.

intensities <- rgamma(1000, 2, 2)
out <- plotCellIntensity(x, y, intensities)
expect_equal(min(intensities), as.double(head(names(out), 1)))
expect_equal(max(intensities), as.double(tail(names(out), 1)))

out <- plotCellIntensity(x, y, intensities, irange=c(1, 4))
expect_equal(1, as.double(head(names(out), 1)))
expect_equal(4, as.double(tail(names(out), 1)))

# Finally, checking we're computing the ranges correctly.

stuff <- list(A=matrix(rgamma(10000, 2, 2), ncol=20))
colnames(stuff$A) <- paste0("X", seq_len(ncol(stuff$A)))
cd <- prepareCellData(stuff)

x <- intensityRanges(cd, p=0.01)
ref <- apply(stuff$A, 2, quantile, p=c(0.01, 0.99))
rownames(ref) <- c("min", "max")
expect_equal(x, ref)

x <- intensityRanges(cd, p=0.05)
ref <- apply(stuff$A, 2, quantile, p=c(0.05, 0.95))
rownames(ref) <- c("min", "max")
expect_equal(x, ref)

#####################################################
# Testing the function that finds my neighbours.
# Couldn't be bothered putting this in a separate file.

set.seed(400)

for (tol in c(0.3,0.4,0.5,0.6,0.7)) { 
    stuff <- list(A=matrix(rgamma(1000, 2, 2), ncol=20))
    colnames(stuff$A) <- paste0("X", seq_len(ncol(stuff$A)))
    cd <- prepareCellData(stuff)
    restuff <- matrix(rgamma(10000, 2, 2), nrow=20)

    # Standard calculation.
    distance <- tol*sqrt(nrow(restuff)) 
    ci <- cellIntensities(cd)
    was.counted <- .Call(cydar:::cxx_find_counted, ci, metadata(cd)$cluster.centers, metadata(cd)$cluster.info, restuff, distance)

    # Reference calculation.
    is.okay <- logical(ncol(restuff))
    for (i in seq_len(ncol(restuff))){ 
        all.distances <- sqrt(colSums((restuff[,i] - ci)^2))
        is.okay[i] <- any(all.distances <= distance)
    }
    expect_identical(is.okay, was.counted)
}


