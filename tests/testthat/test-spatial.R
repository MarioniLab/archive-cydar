#####################################################
# Testing functions related to computing the spatial FDR.

require(cydar); require(testthat)
set.seed(300)
ncells <- 1000
nmarkers <- 10
coords <- matrix(rnorm(ncells*nmarkers, sd=1), nrow=ncells, ncol=nmarkers)
colnames(coords) <- paste0("X", seq_len(nmarkers))
suppressWarnings(cd <- prepareCellData(list(A=coords)))

for (nn in c(10L, 20L, 50L)) { 
    # Testing the nearest neighbour machinery.
    ci <- cellIntensities(cd)
    allbands <- .Call(cydar:::cxx_get_knn_distance, ci, metadata(cd)$cluster.centers, metadata(cd)$cluster.info, nn)
    
    refdist <- as.matrix(dist(t(ci)))
    refbands <- apply(refdist, 1, function(x) { sort(x)[nn] })
    names(refbands) <- NULL
    expect_equal(allbands, refbands)

    # This tests the density calculation machinery.
    bandwidth <- median(refbands)
    densities <- .Call(cydar:::cxx_compute_density, ci, metadata(cd)$cluster.centers, metadata(cd)$cluster.info, bandwidth)

    weightmat <- 1 - (refdist/bandwidth)^3
    weightmat[weightmat < 0] <- 0
    refdens <- rowSums(weightmat^3)
    names(refdens) <- NULL
    expect_equal(refdens, densities)

    # Testing the overall behaviour.
    pval <- rbeta(ncells, 1, 10)
    x <- spatialFDR(coords, pval, neighbors=nn)
    x2 <- spatialFDR(coords, pval, bandwidth=bandwidth)
    expect_equal(x, x2)
    x3 <- spatialFDR(coords, pval, neighbors=nn, naive=TRUE)
    expect_equal(x, x2)

    # Comparing to a straight-up implementation.
    refdist <- as.matrix(dist(coords))
    refbands <- apply(refdist, 1, function(x) { sort(x)[nn] })
    weightmat <- 1 - (refdist/bandwidth)^3
    weightmat[weightmat < 0] <- 0
    refdens <- rowSums(weightmat^3)

    o <- order(pval)
    pval <- pval[o]
    w <- 1/refdens
    w <- w[o]
    qval <- rev(cummin(rev(pval * sum(w)/cumsum(w))))
    qval[o] <- pmin(1, qval)
    names(qval) <- NULL
    expect_equal(x, qval)

    # Checking what happens when we add NAs in the columns.
    discarded <- c(1,4,5,6,7)
    subcoords <- coords
    subcoords[,discarded] <- NA
    subx <- spatialFDR(subcoords, pval, neighbors=nn)
    subx.2 <- spatialFDR(coords[,-discarded], pval, neighbors=nn)
    expect_equal(subx, subx.2)

    # Checking what happens when we add NAs in the pvalues.
    pval.2 <- pval
    discarded <- c(10, 15, 20, 50)
    pval.2[discarded] <- NA_real_
    nax <- spatialFDR(coords, pval.2, neighbors=nn)
    nax.2 <- spatialFDR(coords[-discarded,], pval[-discarded], neighbors=nn)
    expect_equal(nax[-discarded], nax.2)
    expect_true(all(is.na(nax[discarded])))
}

# Testing silly outputs.
pval <- rbeta(ncells, 1, 10)
suppressWarnings(x <- spatialFDR(coords, pval, neighbors=1))
suppressWarnings(x2 <- spatialFDR(coords, pval, bandwidth=0))
expect_equal(x, x2)
expect_equal(x, p.adjust(pval, method="BH"))

