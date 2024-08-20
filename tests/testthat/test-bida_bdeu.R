test_that("posterior_mean.bida_bdeu", {
  nlev <- 2:5
  data <- sapply(nlev, sample, size = 100, replace = T) -1
  y <- 1
  x <- 2
  z <- 3:4
  ess <- 1

  # expected outcome
  nyxz <- counts_from_data_matrix(data, nlev, FALSE) + ess/prod(nlev)
  exp  <- nyxz/rep(colSums(nyxz), each = dim(nyxz)[1])

  bdeu <- bida_bdeu(data, y, c(x, z), 1, nlev)
  obj  <- posterior_mean(bdeu, nlev[x])
  expect_equal(obj, exp, tolerance = 10e-10)
})


test_that("backdoor_mean.bida_bdeu", {

  nlev <- 2:5
  data <- sapply(nlev, sample, size = 10, replace = T) -1
  y <- 1
  x <- 2
  z <- 3:4
  ess <- 1

  # expected outcome
  ayxz <- counts_from_data_matrix(data, nlev, FALSE) + ess/prod(nlev)
  dim(ayxz) <- c(nlev[c(y, x)], prod(nlev[z]))
  axz  <- colSums(ayxz)
  az   <- colSums(axz)

  px.z  <- axz/rep(az, each = nlev[x])
  exp  <- rowSums(sweep(ayxz, 2:3, px.z, "/"), dims = 2)/sum(az)

  bdeu <- bida_bdeu(data, y, c(x, z), 1, nlev)
  obj  <- backdoor_mean(bdeu, nlev[x])
  expect_equal(obj, exp, tolerance = 10e-10)

  obj  <- rowMeans(backdoor_sample(bdeu, 10**4), dims = 2)
  expect_equal(obj, exp, tolerance = 10e-3)
})
