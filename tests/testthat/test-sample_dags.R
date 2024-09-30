test_that("multiplication works", {

  # create bn where one node has many parents
  n <- 3
  dag <- matrix(0, n, n)
  dag[-n, n] <- 1
  nlev = rep(2, n)
  ess <- 1

  set.seed(1)
  bn <- rand_bn(dag, "cat", nlev = nlev, alpha = 1)
  data <- sample_data_from_bn(bn, 10000)
  scorepar <- define_scoreparameters(data, "bdecat", par = list(nlev = nlev, ess = ess))

  #  init search space retrns full DAG
  hardlimit  <- 1
  startspace <-  init_search_space(scorepar, "pcskel", hardlimit = hardlimit, verbose = F)
  expect_true(all(colSums(startspace) <= hardlimit))
})
