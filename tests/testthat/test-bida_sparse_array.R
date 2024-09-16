test_that("multiplication works", {

  x <- bida_sparse_array(0:8, 0:8, 9, default = 0)
  expect_equal(as.array(x), array(0:8, 9))

  # dims ----
  dims   <- c(3, 3, 3)
  names <- list(x = 0:2, y = 0:2, z = 0:2)
  dim(x) <- dims
  dimnames(x) <- names
  expect_equal(dim(x), dims)
  expect_equal(dimnames(x), names)

  # sums ----
  expect_equal(colSums(as.array(x)), as.array(colSums(x)))
  expect_equal(rowSums(as.array(x)), c(as.array(rowSums(x)))) # base::*Sums remove dimnames

  x$default <- 1
  expect_equal(colSums(as.array(x), dims = 2), c(as.array(colSums(x, dims = 2))))
  expect_equal(rowSums(as.array(x), dims = 2), as.array(rowSums(x, dims = 2)))

  # arithmetics
  y <- bida_sparse_array(0:2, 0:2, dim = dim(x), default = 2)
  expect_equal(as.array(x+y), as.array(x)+as.array(y))
  expect_equal(as.array(x-y), as.array(x)-as.array(y))
  expect_equal(as.array(x*y), as.array(x)*as.array(y))
  expect_equal(as.array(x/y), as.array(x)/as.array(y))

})
