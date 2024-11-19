test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("bida_posterior_local_bdeu returns correct counts", {

  n <- 3
  nlev <- rep(3, n)
  data <- expand_grid_fast(k = nlev)

  ps  <- matrix(list(), n)
  ps[[2]] <- list(list(1L), 1)
  bida_posteriors <- bida_posterior_local_bdeu(ps, data, 2, 1:n, 1, nlev)

  # effect of 2 on 1
  az <- array(9, 3) + 1/3
  expected <- new_bida_posterior_bdeu(list(az), 1, 1, c(3, 3))
  expect_equal(bida_posteriors[[1]], expected)

  # effect of 2 on 2
  expect_equal(bida_posteriors[[2]], NULL)

  # effect of 2 on 3
  indx <- data[, 3] + 3*data[, 2] + 3**2*data[, 1]
  ayxz <- bida_sparse_array(rep(1, length(indx)), indx, nlev) + 1/prod(nlev)
  expected <- new_bida_posterior_bdeu(list(ayxz), 1, 0, c(3, 3))
  expect_equal(bida_posteriors[[3]], expected)
})

test_that("methods for computing means and sample", {
  n <- 3
  nlev <- rep(3, n)
  data <- expand_grid_fast(k = nlev)

  ps  <- list(list(1L), 1)
  bida_posteriors <- bida_posterior_local_bdeu(ps, data, 2, 1:n, 1, nlev)

  # check that means are correct
  mean <- array(1/3, dim = c(3, 3))
  expect_equal(mean.bida_posterior_bdeu(bida_posteriors[[1]]), mean)
  expect_equal(mean.bida_posterior_bdeu(bida_posteriors[[3]]), mean)

  # check sample has correct dimensions
  N <- 10
  dims <- c(3, 3, N)
  smpl <- sample.bida_posterior_bdeu(bida_posteriors[[1]], N)
  expect_equal(dim(smpl), dims)

  smpl <- sample.bida_posterior_bdeu(bida_posteriors[[3]], N)
  expect_equal(dim(smpl), dims)

  smpl <- sample.bida_posterior_bdeu(bida_posteriors[[1]], N, contrasts = list(jsd = JSD))
  expect_equal(dim(smpl), c(N, 1))
  expect_equal(smpl, matrix(0, nrow = N, dimnames = list(NULL, c("jsd"))))

  smpl <- sample.bida_posterior_bdeu(bida_posteriors[[3]], N, contrasts = list(jsd = JSD))
  expect_equal(dim(smpl), c(N, 1))

})

