test_that("methods works", {

  # specify counts for wich only first level of z is observed
  nlev <- c(2, 2, 3)
  Nyxz <- bida_sparse_array(1:4, 0:3, nlev)
  size <- 10
  ess  <- 1

  # posterior sample
  sample <- posterior_sample.bdeu_posterior(Nyxz, size, ess)
  expect_equal(colSums(sample), array(1, dim = c(nlev[-1], size)))

  sample <- posterior_sample.bdeu_posterior(as.array(Nyxz), size, ess)
  expect_equal(colSums(sample), array(1, dim = c(nlev[-1], size)))

  # backdoor mean
  expect_equal(backdoor_mean.bdeu_posterior(Nyxz, ess),
               backdoor_mean.bdeu_posterior(as.array(Nyxz), ess))

  # backdoor sample
  sample <- backdoor_sample.bdeu_posterior(Nyxz, size, ess)
  expect_equal(colSums(sample), array(1, dim = c(nlev[-c(1, 3)], size)))

  sample <- backdoor_sample.bdeu_posterior(as.array(Nyxz), size, ess)
  expect_equal(colSums(sample), array(1, dim = c(nlev[-c(1, 3)], size)))
})
