test_that("define_scoreparameters", {

  data <- cbind(z = rep(0:2, 2),
                x = rep(c(0, 2), each = 3),
                y = c(rep(0, 3), rep(1, 3)))
  nlev <- rep(3, 3, 3)
  n <- length(nlev)

  df <- data.frame(lapply(seq_len(n),
                          function(x) factor(data[, x], seq.int(0, nlev[x]-1))))
  colnames(df) <- colnames(data)

  # compute counts for sc
  j <- 3
  parentnodes <- seq_along(nlev)[-j]
  counts <- counts_from_data_matrix(data, nlev)
  dim(counts) <- c(prod(nlev[parentnodes]), nlev[j])

  par <- list(ess = 1, edgepf = 2, nlev = nlev)
  pG  <- -length(parentnodes)*log(par$edgepf)
  lookup <- rlang::new_environment()

  # BiDAG-bdecat-score. Note: ignore missing levels
  scorepar <- define_scoreparameters(data, "bdecat", par)
  score <- unname(BiDAG:::DAGcorescore(j, parentnodes, n, scorepar))-pG
  expect_equal(score,
               famscore_bdeu(counts[rowSums(counts)>0, colSums(counts)>0], par$ess))


  # Local-structure
  for (method in c("pcart", "tree", "ptree", "ldag")) {
    par$local_struct = method
    scorepar <- define_scoreparameters(data, "bdecat", par, lookup)
    score <- unname(BiDAG:::usrDAGcorescore(j, parentnodes, n, scorepar))-pG
    tmp <- optimize_partition_from_data(data, 3, 1:2, 1, nlev, par$local_struct, verbose = FALSE)
    expected <-  sum(tmp$scores)
    expect_equal(score, expected)

    # check that bdeu-object is stored in lookup-table
    expect_equal(class(lookup[[par$local_struct]][["3.1.2"]]), "bida_bdeu")

    # check that also second call to score function returns correct score
    # - now this should be returned from the lookup-table
    score <- unname(BiDAG:::usrDAGcorescore(j, parentnodes, n, scorepar))-pG
    expect_equal(score, expected)
  }
})
