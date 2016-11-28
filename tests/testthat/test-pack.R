require(cydar); require(testthat);

# Testing the packing code.
incoming <- as.integer(c(1, 3:5, 7, 8:9, 10))
incoming <- incoming[sample(length(incoming))]
expected <- c(1L, 3L, -5L, 7L, -10L)
expect_identical(packIndices(list(incoming)), list(expected))

incoming2 <- as.integer(c(4:10, 16:20))
expected2 <- c(4L, -10L, 16L, -20L)
expect_identical(packIndices(list(incoming, incoming2)), list(expected, expected2))

incoming3 <- as.integer(c(4:10, 8:20))
expected3 <- c(4L, -20L)
expect_identical(packIndices(list(incoming, incoming2, incoming3)), list(expected, expected2, expected3))

expect_identical(packIndices(list(incoming, integer(0))), list(expected, integer(0)))
expect_identical(packIndices(list()), list())
expect_error(packIndices(list(-1L)), "all indices must be non-negative integers")

# Testing the unpacking code.
incoming <- as.integer(c(3, -5, 10, 14, -18))
expected <- c(3:5, 10L, 14:18)
expect_identical(unpackIndices(list(incoming)), list(expected))

incoming2 <- as.integer(c(10, 18, -20, 30, 40, -45, 60))
expected2 <- c(10L, 18:20, 30L, 40:45, 60L)
expect_identical(unpackIndices(list(incoming, incoming2)), list(expected, expected2))

expect_identical(unpackIndices(list(incoming, integer(0))), list(expected, integer(0)))
expect_identical(unpackIndices(list()), list())

expect_error(unpackIndices(list(-5L)), "inappropriate negative values in compressed index vector")
expect_error(unpackIndices(list(c(1L, -5L, -6L))), "inappropriate negative values in compressed index vector")
expect_error(unpackIndices(list(c(1L, 5L, -3L))), "inappropriate negative values in compressed index vector")
expect_error(unpackIndices(list(c(1L, 5L, 3L))), "absolute values of compressed indices must always increase")
expect_error(unpackIndices(list(0L)), "zero values")

