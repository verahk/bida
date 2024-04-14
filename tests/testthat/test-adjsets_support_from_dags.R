test_that("function return same output as lapply", {

  set.seed(007)
  n <- 10
  adjsets <- c("o_min", "pa_min", "o")
  dags <- replicate(10, randDAG(n, 3), simplify = F)
  tmp <- matrix(0, n, n)
  dags[[length(dags)+1]] <- tmp
  dags[[length(dags)+1]] <- upper.tri(tmp)*1
  support <- rep(1/length(dags), length(dags))

  # compute support for each sample
  system.time(ps <- adjsets_support_from_dags(adjsets, dags))

  #profvis::profvis(adjsets_support_from_dags(adjsets, dags, support))

  # find unique sets using lapply
  system.time(all_sets <- lapply(dags, function(g) adjsets_from_dag(adjsets, g)))

  for (x in seq_len(n)){
    for (y in seq_len(n)[-x]) {
      for (a in adjsets)

      sets <- ps[[a]]$sets[[x, y]]
      supp <- ps[[a]]$support[[x, y]]

      tmp <- lapply(all_sets, function(m) m[[x, y, a]])
      u <- unique(tmp)
      w <- rowsum_fast(support, tmp, u, FALSE)

      # check that each set is sets has a match in u
      for (r in 1:nrow(sets)) {
        z <- sets[r, ]
        z <- z[!is.na(z)]

        pos <- match_vec(z, u)
        expect_false(is.na(pos))
        expect_equal(w[pos], supp[r])
      }
    }
  }
})


test_that("exclude large sets", {

  set.seed(007)
  adjsets <- c("o_min", "pa_min", "o")
  dag <-  upper.tri(diag(4))*1
  dags <- list(dag)
  support <- rep(1/length(dags), length(dags))

  # compute support for each sample
  ps <- adjsets_support_from_dags(adjsets, dags, support, checksize = function(x, y, z)  FALSE)
  expect_equal(ps$o, ps$o_min)
  expect_equal(ps$pa_min, ps$o_min)

  sets <- ps$pa_min$sets
  expect_equal(unique(sets[upper.tri(sets)]), list(matrix(NA)))
  expect_equal(unlist(sets[lower.tri(sets)]), .col(dim(dag))[lower.tri(dag)])

  # compute support for each sample
  ps <- adjsets_support_from_dags(adjsets, dags, support, checksize = function(x, y, z) length(z) <= 1)
  expect_equal(ps$o, ps$o_min)
  expect_equal(ps$o$sets[-4, 4], list(matrix(NA), matrix(1), matrix(NA)))
})
