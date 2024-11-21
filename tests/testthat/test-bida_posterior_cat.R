test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("bida_posterior_cat returns correct counts", {

  n <- 3
  nlev <- rep(3, n)
  ess  <- 1
  data <- expand_grid_fast(k = nlev)

  ps <- list(list(1L), 1)
  x  <- 2
  y  <- 1:n
  bp <- lapply(y, function(yy) bida_posterior_cat(ps, data, x, yy, ess, nlev))
  bp_local  <- bida_posterior_cat_local(ps, data, x, y, ess, nlev)

  # effect of 2 on 1
  Nz <- array(9, 3)
  expected <- new_bida_posterior_cat(list(Nz), 1, 1, ess, c(3, 3))
  expect_equal(bp[[1]], expected)
  expect_equal(bp_local[[1]], expected)

  # effect of 2 on 2
  expect_equal(bp[[2]], NULL)
  expect_equal(bp_local[[2]], NULL)

  # effect of 2 on 3
  indx <- data[, 3] + 3*data[, 2] + 3**2*data[, 1]
  Nyxz <- bida_sparse_array(rep(1, length(indx)), indx, nlev)
  expected <- new_bida_posterior_cat(list(Nyxz), 1, 0, ess, c(3, 3))
  expect_equal(bp[[3]], expected)
  expect_equal(bp_local[[3]], expected)

})

test_that("methods for posterior means and posterior sample", {
  n <- 3
  nlev <- rep(3, n)
  data <- expand_grid_fast(k = nlev)

  ps  <- list(list(1L), 1)
  bida_posteriors <- bida_posterior_cat_local(ps, data, 2, 1:n, 1, nlev)
  means <- lapply(bida_posteriors, posterior_mean)

  # check that means are correct
  mean <- array(1/3, dim = c(3, 3))
  expect_equal(posterior_mean(bida_posteriors[[1]]), mean)
  expect_equal(posterior_mean(bida_posteriors[[3]]), mean)

  # check sample has correct dimensions
  N <- 10
  dims <- c(3, 3, N)
  smpl <- posterior_sample(bida_posteriors[[1]], N)
  expect_equal(dim(smpl), dims)

  smpl <- posterior_sample(bida_posteriors[[3]], N)
  expect_equal(dim(smpl), dims)

  smpl <- posterior_sample(bida_posteriors[[1]], N, contrasts = list(jsd = jsd))
  expect_equal(dim(smpl), c(N, 1))
  expect_equal(smpl, matrix(0, nrow = N, dimnames = list(NULL, c("jsd"))))

  smpl <- posterior_sample(bida_posteriors[[3]], N, contrasts = list(jsd = jsd))
  expect_equal(dim(smpl), c(N, 1))
})

