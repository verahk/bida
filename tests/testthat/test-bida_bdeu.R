test_that("posterior_mean.bida_bdeu", {
  nlev <- 2:5
  data <- sapply(nlev, sample, size = 100, replace = T) -1
  y <- 1
  x <- 2
  z <- 3:4
  ess <- 1

  # expected outcome
  nyxz <- counts_from_data_matrix(data, nlev, FALSE)
  ayxz <- nyxz + ess/prod(nlev)
  exp  <- ayxz/rep(colSums(ayxz), each = dim(nyxz)[1])

  bdeu <- bida_bdeu(data, y, c(x, z), 1, nlev)
  expect_equal(as.array(bdeu$counts), nyxz)
  expect_equal(as.array(colSums(bdeu$counts)), colSums(nyxz))
  expect_equal(posterior_mean(bdeu), exp, tolerance = 10e-10)

  # no parents
  ny <- rowSums(nyxz)
  ay <- ny + ess/nlev[1]
  py <- ay/sum(ay)
  bdeu <- bida_bdeu(data, y, integer(0), 1, nlev)
  expect_equal(as.array(bdeu$counts), ny, ignore_attr = TRUE)
  expect_equal(posterior_mean(bdeu), py, tolerance = 10e-10, ignore_attr = TRUE)
})

test_that("posterior_mean.bida_bdeu with partition", {
  data <- cbind(
    y = c(rep(0, 9), rep(1, 2*9)),
    x = rep(0:2, each = 9),
    z = rep(0:2, 9))
  nlev <- rep(3, 3)
  names(nlev) <- colnames(data)
  ess <- 1
  opt <- optimize_partition_from_data(data, 1, 2:3, nlev, ess = 1, "ptree")

  bdeu <- bida_bdeu(data, 1, 2:3, ess, nlev, partition = opt)
  bdeu2 <- bida_bdeu(data, 1, 2, nlev, ess = 1) # should give same res partition is a binary split
  expect_equal(posterior_mean(bdeu), posterior_mean(bdeu2), ignore_attr = TRUE)
  expect_equal(posterior_mean(bdeu, reduced = FALSE), array(rep(posterior_mean(bdeu2), nlev[3]), nlev))
  expect_equal(backdoor_mean(bdeu), backdoor_mean(bdeu2), ignore_attr = TRUE)
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


