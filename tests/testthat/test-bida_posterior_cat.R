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
  bp <- lapply(y, function(yy) bida_posterior_cat(ps, data, x, yy, ess, nlev)[[yy]])
  bp_local <- bida_posterior_cat(ps, data, x, y, ess, nlev)
  expect_equal(bp, bp_local)
  #bp_local  <- bida_posterior_cat_local(ps, data, x, y, ess, nlev)

  # effect of 2 on 1
  Nz <- array(9, 3)
  expected <- new_bida_posterior_cat(list(Nz), 1, 1, ess, c(3, 3))
  expect_equal(bp[[1]], expected)
  #expect_equal(bp_local[[1]], expected)

  # effect of 2 on 2
  expect_equal(bp[[2]], NULL)
  #expect_equal(bp_local[[2]], NULL)

  # effect of 2 on 3
  indx <- data[, 3] + 3*data[, 2] + 3**2*data[, 1]
  Nyxz <- bida_sparse_array(rep(1, length(indx)), indx, nlev)
  expected <- new_bida_posterior_cat(list(Nyxz), 1, 0, ess, c(3, 3))
  expect_equal(bp[[3]], expected)
  #expect_equal(bp_local[[3]], expected)

})

test_that("posterior_sample with support on zero", {
  n <- 3
  nlev <- 2:4

  dag  <- matrix(0, n, n)
  dag[upper.tri(dag)] <- 1
  dags <- list(dag, t(dag))
  x <- 2
  y <- 3

  data <- sapply(nlev, sample.int, size = 10, replace = TRUE)-1
  ps   <- adjset_support(dag_support(dags), x, y, "o")
  fit  <- bida_posterior_cat(ps, data, 1, 2, nlev = nlev, ess = 1)[[2]]

  smpl <- posterior_sample(fit, 10, contrasts = list(jsd = jsd))
  expect_equal(dim(smpl), c(10, 1))
})

