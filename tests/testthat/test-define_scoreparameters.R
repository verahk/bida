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

  scorepar <- define_scoreparameters(data, "bdecat", par)
  score <- unname(BiDAG:::DAGcorescore(j, parentnodes, n, scorepar))-pG
  expect_equal(score,
               famscore_bdeu(counts[rowSums(counts)>0, colSums(counts)>0], par$ess))


  # Local-structure
  for (method in c("pcart", "tree", "ptree", "treereg")) { # ldag"
    # define score-parameter and compute score
    par$local_struct <- method
    scorepar <- define_scoreparameters(data, "bdecat", par, lookup)
    score <- unname(BiDAG:::usrDAGcorescore(j, parentnodes, n, scorepar))-pG
    opt <- optimize_partition_from_data(data, j, parentnodes, 1, nlev, par$local_struct, verbose = FALSE)
    expect_equal(score, attr(opt, "score"))

    # check that bdeu-object is stored in lookup-table
    class <- ifelse(method == "pcart", method, "tree")
    expect_true(inherits(lookup[[par$local_struct]][["3.1.2"]], class))

    # check score for parent set of size < 2
    score <- unname(BiDAG:::DAGcorescore(j, parentnodes[1], n, scorepar))+log(par$edgepf)
    nodes <- c(parentnodes[1], j)
    tab <- counts_from_data_matrix(data[, nodes], nlev[nodes])
    expect_equal(score, famscore_bdeu(tab, 1))

    score <- unname(BiDAG:::DAGcorescore(j, integer(0), n, scorepar))
    tab <-  matrix(tabulate(data[, j]+1, nlev[j]), ncol = nlev[j])
    expect_equal(score, famscore_bdeu(tab, 1))
  }

})
