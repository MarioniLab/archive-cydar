# This tests the distance calculating machinery.

require(cydar); require(testthat)
set.seed(600)
ncells <- 1000

for (nmarkers in c(5, 10, 20)) { 
    coords <- matrix(rnorm(ncells*nmarkers, sd=1), nrow=ncells, ncol=nmarkers)
    colnames(coords) <- paste0("X", seq_len(nmarkers))
    suppressWarnings(cd <- prepareCellData(list(A=coords)))
    sub.chosen <- c(2,3,4)
    suppressWarnings(cd.sub.1 <- prepareCellData(list(A=coords[,sub.chosen])))
    o.1 <- order(cellData(cd.sub.1)$cell.id)
    suppressWarnings(cd.sub.2 <- prepareCellData(list(A=coords), markers=sub.chosen))
    o.2 <- order(cellData(cd.sub.2)$cell.id)

    for (nn in c(10, 20, 50)) { 
        ci <- cellIntensities(cd)
        refdist <- as.matrix(dist(t(ci)))
        refbands <- apply(refdist, 1, function(x) { sort(x)[2:(nn+1)] })
        dimnames(refbands) <- NULL
        refbands <- t(refbands)

        stuff <- neighborDistances(cd, downsample=1, neighbors=nn, as.tol=FALSE)
        expect_equal(stuff, refbands)
        
        stuff2 <- neighborDistances(cd, downsample=1, neighbors=nn)
        expect_equal(stuff, stuff2*sqrt(nmarkers))

        # Increasing the downsampling.
        ds <- 10
        stuff3 <- neighborDistances(cd, downsample=ds, neighbors=nn)
        chosen <- seq(from=1, to=ncells, by=ds)
        expect_equal(stuff2[chosen,,drop=FALSE], stuff3)

        stuff4 <-  neighborDistances(cd, downsample=ds, neighbors=nn, naive=TRUE)
        expect_equal(stuff3, stuff4)

        # Testing behaviour upon using a subset of markers.
        stuff.sub.1 <- neighborDistances(cd.sub.1, downsample=1, neighbors=nn)
        stuff.sub.2 <- neighborDistances(cd.sub.2, downsample=1, neighbors=nn)
        expect_equal(stuff.sub.1[o.1,], stuff.sub.2[o.2,])
    }
}
