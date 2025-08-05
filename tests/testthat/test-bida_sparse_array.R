test_that("methods applied to standard and sparse arrays give same results", {

  x <- bida_sparse_array(1:9, 0:8, 9, default = 0)
  expect_equal(as.array(x), array(1:9, 9))

  # change dims ----
  dims   <- c(3, 3, 3)
  names <- list(x = 0:2, y = 0:2, z = 0:2)
  dim(x) <- dims
  dimnames(x) <- names
  expect_equal(dim(x), dims)
  expect_equal(dimnames(x), names)

  # as.sparse_array
  y <- bida_sparse_array(as.array(x))
  dimnames(y) <- lapply(dimnames(y), as.integer)  # array() convert all dimnames to character vectors..
  expect_equal(y, x)

  # col- and rowSums ----
  expect_equal(colSums(as.array(x)), as.array(colSums(x)))
  expect_equal(rowSums(as.array(x)), c(as.array(rowSums(x)))) # base::*Sums remove dimnames

  x$default <- 1
  expect_equal(colSums(as.array(x), dims = 2), c(as.array(colSums(x, dims = 2))))
  expect_equal(rowSums(as.array(x), dims = 2), as.array(rowSums(x, dims = 2)))


  # arithmetics
  y <- bida_sparse_array(1:10, 0:9, dim = dim(x), default = 2)
  expect_equal(as.array(x+y), as.array(x)+as.array(y))
  expect_equal(as.array(x-y), as.array(x)-as.array(y))
  expect_equal(as.array(x*y), as.array(x)*as.array(y))
  expect_equal(as.array(x/y), as.array(x)/as.array(y))


  # rep
  expect_equal(c(as.array(rep(x, each = 2))), rep(as.array(x), each = 2))
  expect_equal(c(as.array(rep(x, times = 2))), rep(as.array(x), times = 2))
  expect_equal(c(as.array(rep(x, times = 2, each = 2))), rep(as.array(x), times = 2, each = 2))

  # aperm
  expect_equal(as.array(aperm(x, 3:1)), aperm(as.array(x), 3:1))

  # random index
  xx <- bida_sparse_array(x$value, sample(x$index), c(3, 3))
  yy <- rep(colSums(xx), each = dim(xx)[1])
  expect_equal(as.array(xx), array(replace(rep(0, 9), xx$index+1, xx$value), dim(xx)))
  expect_equal(as.array(yy), rep(colSums(as.array(xx)), each = dim(xx)[1]), ignore_attr = TRUE)
})

test_that("aperm works", {

  # index of length 1
  x <- bida_sparse_array(10L, 3, dim = 1:3, default = 0)
  expect_equal(as.array(aperm(x, c(1, 3, 2))), aperm(as.array(x), c(1, 3, 2)))

  # random index
  x <- bida_sparse_array(1:prod(dim(x)), sample.int(prod(dim(x)))-1, dim = dim(x), default = 0)
  expect_equal(as.array(aperm(x, c(1, 3, 2))), aperm(as.array(x), c(1, 3, 2)))
})
