test_that("multiplication works", {

  # create bn where one node has many parents
  n <- 3
  dag <- matrix(0, n, n)
  dag[-n, n] <- 1
  nlev = rep(2, n)
  ess <- 1

  set.seed(1)
  bn <- rand_bn(n, 2, "cat", nlev = rep(2, n))
  data <- sample_data_from_bn(bn, 10000)
  scorepar <- define_scoreparameters(data, "bdecat", par = list(nlev = nlev, ess = ess))

  #  init search space retrns full DAG
  hardlimit  <- 1
  expect_warning(init_search_space(scorepar, "pcskel", hardlimit = hardlimit, verbose = T))

  if (FALSE) {
    smpl <- sample_dags(scorepar, hardlimit = hardlimit, verbose = T)
    npar <- sapply(lapply(smpl$traceadd$incidence, as.matrix), colSums)
    expect_true(max(npar) <= hardlimit+1) # partitionMCMC can add 1 parent with respect to startspace
  }
})
