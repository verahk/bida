
test_that("optimization routines handles mixed cardinality data", {


  levels <- list(0:1, 0:2, 0:3)
  r <- 3
  q <- prod(lengths(levels))


  # construct counts that collapses all rows of the CPT ----
  counts <- cbind(1, 10, rep(10**2, q))


  for (m in c("tree", "ldag")) {
    res <- optimize_partition(counts, levels, ess = 1, method = m, regular = F, verbose = F)
    expect_true(length(res$partition) == 1)
  }


  res <- optimize_partition(counts, levels, ess = 1, method = "tree", regular = T, verbose = F)
  expect_true(length(res$partition) == q)

  res <- optimize_partition(counts, levels, ess = 1, method = "ldag", regular = T, verbose = F)
  #cbind(expand.grid(levels), counts, get_parts(res$partition))
  expect_true(length(res$partition) == 2)


  # construct 2-region split counts ----
  counts <- cbind(1, 10, c(rep(10, q-1), 1000))

  res <- optimize_partition(counts, levels, ess = 1, method = "tree", regular = F, verbose = F)
  expect_true(length(res$partition) == 7)

  res <- optimize_partition(counts, levels, ess = 1, method = "tree", regular = T, verbose = F)
  expect_true(length(res$partition) == 7)

  res <- optimize_partition(counts, levels, ess = 1, method = "ldag", regular = F, verbose = F)
  expect_true(length(res$partition) == 2)

  res <- optimize_partition(counts, levels, ess = 1, method = "ldag", regular = T, verbose = F)
  expect_true(length(res$partition) == 2)
})



test_that("regularity options work as intended and gives correct scores", {

  # counts that returns independence
  levels <- list(0:1, 0:1)
  counts <- cbind(1, rep(10, 4))

  # regular = FALSE ----
  regular <- FALSE

  ## expected values
  part_exp  <- rep(1, 4)
  score_exp <- famscore_bdeu_1row(colSums(counts), ess = 1, r = 2, q = 1, s = 1)

  ### tree
  res <- optimize_partition(counts, levels, ess = 1, method = "tree", regular = regular)
  expect_equal(get_parts(res$partition), part_exp, ignore_attr = T)
  expect_equal(res$scores, score_exp, ignore_attr = T)


  ### ldag
  res <- optimize_partition(counts, levels, ess = 1, method = "ldag", regular = regular)
  expect_equal(get_parts(res$partition), part_exp, ignore_attr = T)
  expect_equal(res$scores, score_exp, ignore_attr = T)

  ### part
  res <- optimize_partition(counts, levels, ess = 1, method = "part", regular = regular)
  expect_equal(get_parts(res$partition), part_exp, ignore_attr = T)
  expect_equal(res$scores, score_exp, ignore_attr = T)


  # regular = TRUE ----
  regular <- TRUE

  ## expected values - for tree
  part_exp  <- 1:4
  score_exp <- sum(famscore_bdeu_byrow(counts, ess = 1, r = 2, q = 4, s = 1))

  ### tree
  res <- optimize_partition(counts, levels, ess = 1, method = "tree", regular = regular)
  expect_true(all(get_parts(res$partition) %in% part_exp))
  expect_equal(sum(res$scores), score_exp, ignore_attr = T)


  ## expected values - for ldag and part
  part_exp  <- 1:2
  score_exp <- sum(famscore_bdeu_1row(colSums(counts[-1, ]), ess = 1, r = 2, q = 4, s = 3))
  score_exp <- score_exp + famscore_bdeu_1row(counts[1, ], ess = 1, r = 2, q = 4, s = 1)

  ### ldag
  res <- optimize_partition(counts, levels, ess = 1, method = "ldag", regular = regular)
  expect_true(all(get_parts(res$partition) %in% part_exp))
  expect_equal(sum(res$scores), score_exp, ignore_attr = T)

  ### part
  res <- optimize_partition(counts, levels, ess = 1, method = "part", regular = regular)
  expect_true(all(get_parts(res$partition) %in% part_exp))
  expect_equal(sum(res$scores), score_exp, ignore_attr = T)


})
