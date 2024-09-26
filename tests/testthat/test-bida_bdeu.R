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

test_that("aperm.bida_bdeu works", {

  data <- cbind(rep(0, 10), rep(1, 10), rep(2, 10))
  nlev <- c(2, 3, 4)
  col <- seq_len(prod(nlev[-1]))-1
  parts <- rep(seq_len(nlev[3]), each = nlev[2])  # split on parent number 2
  partition <- split(col, parts)

  bdeu <- bida_bdeu(data, 1, 2:3, 1, nlev, partition)
  bdeu_perm <- aperm(bdeu, c(1, 3, 2))
  partition_perm <- split(col%%nlev[2]*nlev[3]+col%/%nlev[2]%%nlev[3], parts)

  expect_equal(as.array(bdeu_perm$counts), aperm(as.array(bdeu$counts), c(1, 3, 2)))
  expect_equal(bdeu_perm$partition, partition_perm)
})


