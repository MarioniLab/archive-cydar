# This tests the construction and function of the cyData class.

require(testthat); require(cydar); require(S4Vectors)

set.seed(1000)
counts <- matrix(rpois(1000, 10), ncol=10)
medians <- matrix(rgamma(2000, 1, 1), ncol=20)
cell.int <- matrix(rgamma(10000, 1, 1), nrow=20)
cell.data <- DataFrame(sample.id=rep(1, 500))
marker.data <- DataFrame(row.names=LETTERS[1:20])

# Checking constructor.
cyd <- cyData(counts, markerData=marker.data, intensities=medians, cellIntensities=cell.int, cellData=cell.data)
expect_error(cyData(counts[1:10,], markerData=marker.data, intensities=medians, cellIntensities=cell.int, cellData=cell.data),
             "'intensities' and 'object' must have same number of rows")
expect_error(cyData(counts, markerData=marker.data[1:10,], intensities=medians, cellIntensities=cell.int, cellData=cell.data),
             "'markerData' and 'cellIntensities' must have same number of rows")
expect_error(cyData(counts, markerData=marker.data, intensities=medians[1:10,], cellIntensities=cell.int, cellData=cell.data),
             "'intensities' and 'object' must have same number of rows")
expect_error(cyData(counts, markerData=marker.data, intensities=medians[,1:10], cellIntensities=cell.int, cellData=cell.data),
             "number of rows in 'markerData' and columns in 'intensities' must be equal")
expect_error(cyData(counts, markerData=marker.data, intensities=medians, cellIntensities=cell.int[1:10,], cellData=cell.data),
             "'markerData' and 'cellIntensities' must have same number of rows")
expect_error(cyData(counts, markerData=marker.data, intensities=medians, cellIntensities=cell.int, cellData=cell.data[1:10,,drop=FALSE]),
             "number of rows in 'cellData' and columns in 'cellIntensities' should be equal")

# Checking getters.
expect_identical(marker.data, markerData(cyd))
expect_identical(medians, intensities(cyd))
expect_identical(cell.int, cellIntensities(cyd))
expect_identical(cell.data, cellData(cyd))
expect_equivalent(counts, assay(cyd))
expect_identical(nrow(marker.data), nmarkers(cyd))       
expect_identical(ncol(cell.int), ncells(cyd))      

# Checking that arguments are passed to SummarizedExperiment.
cyd <- cyData(counts, metadata=list("YAY"), markerData=marker.data, intensities=medians, cellIntensities=cell.int, cellData=cell.data)
expect_identical("YAY", metadata(cyd)[[1]])

# Running setters.
cyd.x <- cyd
random.values <- runif(nmarkers(cyd))
markerData(cyd.x)$score <- random.values
expect_identical(markerData(cyd.x)$score, random.values)

new.names <- paste0("X", seq_len(nmarkers(cyd)))
rownames(markerData(cyd.x)) <- new.names
expect_identical(new.names, rownames(markerData(cyd.x)))

expect_error(markerData(cyd.x) <- marker.data[1:10,], "'markerData' and 'cellIntensities' must have same number of rows")

cell.ids <- seq_len(500)
cellData(cyd.x)$cell.id <- cell.ids
expect_identical(cellData(cyd.x)$cell.id, cell.ids)
expect_error(cellData(cyd.x) <- cell.data[1:10,,drop=FALSE], "number of rows in 'cellData' and columns in 'cellIntensities' should be equal")

more.cell.int <- matrix(rgamma(10000, 1, 1), nrow=20)
cellIntensities(cyd.x) <- more.cell.int
expect_identical(cellIntensities(cyd.x), more.cell.int)
expect_error(cellIntensities(cyd.x) <- cell.int[1:10,], "'markerData' and 'cellIntensities' must have same number of rows")

more.medians <- matrix(rgamma(2000, 1, 1), ncol=20)
intensities(cyd.x) <- more.medians
expect_identical(intensities(cyd.x), more.medians)
expect_error(intensities(cyd.x) <- medians[,1:10], "number of rows in 'markerData' and columns in 'intensities' must be equal")

# Running subsetters.
by.row <- c(10, 20, 50, 60, 71, 92, 100)
cyd.x <- cyd[by.row,]
expect_identical(assay(cyd.x), assay(cyd)[by.row,])
expect_identical(intensities(cyd.x), intensities(cyd)[by.row,])
expect_identical(cellIntensities(cyd.x), cellIntensities(cyd))
expect_identical(cellData(cyd.x), cellData(cyd))
expect_identical(markerData(cyd.x), markerData(cyd))

by.col <- c(1, 2, 5, 9)
cyd.x <- cyd[,by.col]
expect_identical(assay(cyd.x), assay(cyd)[,by.col])
expect_identical(intensities(cyd.x), intensities(cyd))
expect_identical(cellIntensities(cyd.x), cellIntensities(cyd))
expect_identical(cellData(cyd.x), cellData(cyd))
expect_identical(markerData(cyd.x), markerData(cyd))

cyd.x <- cyd[by.row, by.col]
expect_identical(assay(cyd.x), assay(cyd)[by.row,by.col])
expect_identical(intensities(cyd.x), intensities(cyd)[by.row,])
expect_identical(cellIntensities(cyd.x), cellIntensities(cyd))
expect_identical(cellData(cyd.x), cellData(cyd))
expect_identical(markerData(cyd.x), markerData(cyd))

# Running subset replacement.
counts2 <- matrix(rpois(1000, 10), ncol=10)
medians2 <- matrix(rgamma(2000, 1, 1), ncol=20)
cyd2 <- cyData(counts2, markerData=marker.data, intensities=medians2, cellIntensities=cell.int, cellData=cell.data)
cyd3 <- cyData(counts2, markerData=marker.data[rev(seq_len(nrow(marker.data))),], 
    intensities=medians2, cellIntensities=cell.int, cellData=cell.data)

cyd.x <- cyd
cyd.x[by.row,] <- cyd2[by.row,]
counts.x <- counts
counts.x[by.row,] <- counts2[by.row,]
expect_equivalent(counts.x, assay(cyd.x))
medians.x <- medians
medians.x[by.row,] <- medians2[by.row,]
expect_identical(medians.x, intensities(cyd.x))
expect_identical(cellIntensities(cyd.x), cellIntensities(cyd))
expect_identical(cellData(cyd.x), cellData(cyd))
expect_identical(markerData(cyd.x), markerData(cyd))

cyd.x <- cyd
cyd.x[,by.col] <- cyd2[,by.col]
counts.x <- counts
counts.x[,by.col] <- counts2[,by.col]
expect_equivalent(counts.x, assay(cyd.x))
expect_identical(medians, intensities(cyd.x))
expect_identical(cellIntensities(cyd.x), cellIntensities(cyd))
expect_identical(cellData(cyd.x), cellData(cyd))
expect_identical(markerData(cyd.x), markerData(cyd))

expect_error(cyd.x[by.row,] <- cyd3[by.row,], "'markerData' are not identical")
expect_error(cyd.x[,by.col] <- cyd3[,by.col], "'markerData' are not identical")

# Running combiners.
cyd.x <- cbind(cyd, cyd)
expect_equivalent(assay(cyd.x), cbind(assay(cyd), assay(cyd)))
expect_equivalent(markerData(cyd.x), markerData(cyd))
expect_equivalent(intensities(cyd.x), intensities(cyd))
expect_equivalent(cellIntensities(cyd.x), cellIntensities(cyd))
expect_identical(cellData(cyd.x), cellData(cyd))
expect_error(cbind(cyd, cyd2), "'intensities' are not identical")
expect_error(cbind(cyd, cyd3), "'markerData' are not identical")

cyd.x <- rbind(cyd, cyd2)
expect_equivalent(assay(cyd.x), rbind(assay(cyd), assay(cyd2)))
expect_equivalent(intensities(cyd.x), rbind(intensities(cyd), intensities(cyd2)))
expect_equivalent(markerData(cyd.x), markerData(cyd))
expect_equivalent(cellIntensities(cyd.x), cellIntensities(cyd))
expect_identical(cellData(cyd.x), cellData(cyd))
expect_error(rbind(cyd, cyd3), "'markerData' are not identical")

# Checking that empty construction still works.
cyd.minimal <- cyData(assays=matrix(0, 0, 10), markerData=marker.data, cellIntensities=cell.int, cellData=cell.data)
expect_equal(rbind(cyd, cyd.minimal), cyd)
cyd.minimal <- cyData(assays=matrix(0, nrow(medians), 0), intensities=medians, 
    markerData=marker.data, cellIntensities=cell.int, cellData=cell.data)
expect_equal(cbind(cyd, cyd.minimal), cyd)

expect_identical(nrow(cyd[0,]), 0L)
expect_identical(ncol(cyd[,0]), 0L)

