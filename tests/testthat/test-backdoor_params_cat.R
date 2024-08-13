
test_that("constructor works", {
  nlev <- c(2, 3, 3, 3, 3)
  lev  <- lapply(nlev-1, seq.int, from =  0)
  data <- expand_grid_fast(lev)

  # dim of x, y, z
  k <- c(nlev[c(1, 2)], prod(nlev[3:5]))


  # all levels observed
  counts <- array(1, k)
  par <- backdoor_params_cat(data, 2, 1, 3:5, nlev)
  expect_equal(par, counts, ignore_attr = T)
  expect_equal(dim(par), k, ignore_attr = T)


  # some adjustment vars not observed
  data2 <- data[!data[, 5] == 2, ]
  nobs  <- prod(nlev)-prod(nlev[-5])
  counts <- c(rep(1, nobs), rep(0, prod(nlev)-nobs))

  par <- backdoor_params_cat(data2, 2, 1, 3:5, nlev, min_sparse = Inf)
  expect_equal(par, counts, ignore_attr = T)
  expect_equal(dim(par), k, ignore_attr = T)

  par <- backdoor_params_cat(data2, 2, 1, 3:5, nlev, min_sparse = 0)
  expect_equal(par$value, counts[counts > 0], ignore_attr = T)
  expect_equal(par$dim, k, ignore_attr = T)


  # no adjustment vars
  par <- backdoor_params_cat(data, 2, 1, integer(0), nlev)
  counts <- array(k[3], dim = k[1:2])
  expect_equal(par, counts, ignore_attr = T)
  expect_equal(dim(par), k[1:2], ignore_attr = T)

  # adjustment set include y
  par <- backdoor_params_cat(data, 2, 1, c(1, 5), nlev)
  counts <- array(prod(k[-1]), k[1])
  expect_equal(par, counts, ignore_attr = T)
  expect_equal(dim(par), k[1], ignore_attr = T)
})

