# Testing the prepareCellData machinery.

require(testthat); require(cydar)

# Setup.
nmarkers <- 10

ncells1 <- 1001
all.values1 <- matrix(rnorm(ncells1*nmarkers, sd=1), nrow=ncells1, ncol=nmarkers)
colnames(all.values1) <- paste0("X", seq_len(nmarkers))

ncells2 <- 2001
all.values2 <- matrix(rnorm(ncells2*nmarkers, sd=1), nrow=ncells2, ncol=nmarkers)
colnames(all.values2) <- colnames(all.values1)
fs <- list(A=all.values1, B=all.values2)

out <- prepareCellData(list(X=all.values1, Y=all.values2))

# Checking that the intensities are the same.
current <- paste0(cellData(out)$sample.id, ".", cellData(out)$cell.id)
original <- paste0(rep(1:2, c(ncells1, ncells2)), ".", c(seq_len(ncells1), seq_len(ncells2)))
m <- match(original, current)
expect_equivalent(rbind(all.values1, all.values2), t(cellIntensities(out))[m,])

# Checking cluster IDs.
cluster.info <- metadata(out)$cluster.info
all.center.dex <- sapply(cluster.info, "[[", i=1) + 1
all.center.dist <- lapply(cluster.info, "[[", i=2)
all.cluster.centers <- metadata(out)$cluster.centers[,rep(seq_along(all.center.dist), lengths(all.center.dist))]
actual.dist <- sqrt(colSums((all.cluster.centers - cellIntensities(out))^2))
expect_equal(unname(actual.dist), unlist(all.center.dist))
expect_true(all(!sapply(all.center.dist, is.unsorted)))

# Trying out the naive approach.
nout <- prepareCellData(list(X=all.values1, Y=all.values2), naive=TRUE)
expect_equivalent(t(cellIntensities(nout)), rbind(all.values1, all.values2))
expect_identical(cellData(nout)$sample.id, rep(1:2, c(ncells1, ncells2)))
expect_identical(cellData(nout)$cell.id, c(seq_len(ncells1), seq_len(ncells2)))

expect_identical(NULL, metadata(nout)$cluster.info)
expect_identical(NULL, metadata(nout)$cluster.centers)
