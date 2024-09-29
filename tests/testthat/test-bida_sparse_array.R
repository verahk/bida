test_that("multiplication works", {

  x <- bida_sparse_array(1:9, 0:8, 9, default = 0)
  expect_equal(as.array(x), array(0:8, 9))

  # change dims ----
  dims   <- c(3, 3, 3)
  names <- list(x = 0:2, y = 0:2, z = 0:2)
  dim(x) <- dims
  dimnames(x) <- names
  expect_equal(dim(x), dims)
  expect_equal(dimnames(x), names)

  # as.sparse_array
  y <- as.bida_sparse_array(as.array(x))
  dimnames(y) <- lapply(dimnames(y), as.integer)  # array() convert all dimnames to character vectors..
  expect_equal(y, x)


  # arithmetics
  y <- bida_sparse_array(1:2, range(x$index), dim = dim(x), default = 2)
  expect_equal(as.array(x+y), as.array(x)+as.array(y))
  expect_equal(as.array(x-y), as.array(x)-as.array(y))
  expect_equal(as.array(x*y), as.array(x)*as.array(y))
  expect_equal(as.array(x/y), as.array(x)/as.array(y))


  # col- and rowSums ----
  expect_equal(colSums(as.array(x)), as.array(colSums(x)))
  expect_equal(rowSums(as.array(x)), c(as.array(rowSums(x)))) # base::*Sums remove dimnames

  x$default <- 1
  expect_equal(colSums(as.array(x), dims = 2), c(as.array(colSums(x, dims = 2))))
  expect_equal(rowSums(as.array(x), dims = 2), as.array(rowSums(x, dims = 2)))

  # rep
  expect_equal(c(as.array(rep(x, each = 2))), rep(as.array(x), each = 2))
  expect_equal(c(as.array(rep(x, times = 2))), rep(as.array(x), times = 2))
  expect_equal(c(as.array(rep(x, times = 2, each = 2))), rep(as.array(x), times = 2, each = 2))

  # asplit
  y <- as.array(x)
  for (i in seq_len(dim(x)[2])) expect_equal(as.array(asplit(x, 2)[[i]]), asplit(y, 2)[[i]])
  for (i in seq_len(prod(dim(x)[1:2]))) expect_equal(as.array(asplit(x, 1:2)[[i]]), asplit(y, 1:2)[[i]])

  # aperm
  expect_equal(as.array(aperm(x, 3:1)), aperm(as.array(x), 3:1))
})
