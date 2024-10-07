test_that("function works on binary split data with missing variables", {
  data <- cbind(
    y = c(rep(0, 9), rep(1, 2*9)),
    x = rep(0:2, each = 9),
    z = rep(0:2, 9))
  nlev <- rep(3, 3)
  names(nlev) <- colnames(data)

  tab <- matrix(counts_from_data_matrix(data, nlev), byrow = TRUE, ncol = nlev[1])
  ess <- 1
  newdata <- expand.grid(lapply(nlev[-1]-1, seq.int, from = 0))

  # tree
  fit <- optimize_partition_from_data_tree(data, ess, nlev, min_improv = 0, FALSE, FALSE, FALSE)
  parts <- rep(1:3, times = 3)
  sizes <- tabulate(parts)
  scores <- famscore_bdeu_byrow(rowsum(tab, parts), ess, nlev[1], prod(nlev[-1]), sizes)
  expect_equal(predict(fit, newdata), parts)
  expect_equal(attr(fit, "score"), sum(scores))
  expect_equal(attr(fit, "sizes"), sizes)
  expect_true(all(table(predict(fit, newdata) %in% c(0, 3))))
  expect_equal(fit, optimize_partition_from_data(data, 1, 2:3, 1, nlev, "tree"))

  # regular tree
  fit <- optimize_partition_from_data_tree(data, ess, nlev, min_improv = 0, FALSE, TRUE, FALSE)
  parts <- 1:prod(nlev[-1])
  sizes <- tabulate(parts)
  scores <- famscore_bdeu_byrow(rowsum(tab, parts), ess, nlev[1], prod(nlev[-1]), sizes)
  expect_equal(length(unique(predict(fit, newdata))), length(unique(parts)))
  expect_equal(attr(fit, "score"), sum(scores))
  expect_equal(attr(fit, "sizes"), sizes)
  expect_equal(fit, optimize_partition_from_data(data, 1, 2:3, 1, nlev, "treereg"))

  # pcart
  df <- data.frame(lapply(asplit(data, 2), as.factor))
  levels(df$y) <- 0:2
  fit <- optimize_partition_from_df_pcart(df, ess)
  parts <- replace(rep(2, 9), seq(1, 9, by = 3), 1)
  sizes <- tabulate(parts)
  scores <- famscore_bdeu_byrow(rowsum(tab, parts), ess, nlev[1], prod(nlev[-1]), sizes)
  expect_equal(predict(fit, newdata), parts)
  expect_equal(attr(fit, "score"), sum(scores))
  expect_equal(attr(fit, "sizes"), sizes)
  expect_equal(fit, optimize_partition_from_data(data, 1, 2:3, 1, nlev, "pcart"))
})

test_that("pruning procedure and summary procedure works", {
  set.seed(007)
  p <- c(.9, .1, .1, .9)
  z <- runif(1000) > .5
  x <- runif(1000) > .5
  y <- runif(1000) > p[z + 2*x+1]

  data <- cbind(y = y, z = z, x = x)*1
  nlev <- c(2, 2, 2)
  newdata <- expand_grid_fast(k = nlev[-1])
  colnames(newdata) <- colnames(data)[-1]

  tab <- matrix(counts_from_data_matrix(data, nlev), byrow = TRUE, ncol = nlev[1])
  ess <- 1

  # tree: greedy decision tree
  fit <- optimize_partition_from_data(data, 1, 2:3, ess, nlev, "tree", verbose = FALSE)
  summ <- summary(fit)
  expect_equal(summ[, "part"], 1, ignore_attr = T)
  expect_equal(summ[, "score"],
               round(famscore_bdeu_1row(tabulate(data[, 1]+1), 1, 2), 2),
               ignore_attr = T)
  expect_equal(summ[, "size"], 4, ignore_attr = T)


  fit <- optimize_partition_from_data(data, 1, 2:3, ess, nlev, "ptree", verbose = F)
  summ <- summary(fit)
  expect_equal(summ[, "part"], 1:4, ignore_attr = T)
  parts <- predict(fit, newdata)
  expect_true(all(famscore_bdeu_byrow(rowsum(tab, parts), 1, 2, 4, tabulate(parts)) %in% famscore_bdeu_byrow(tab, 1)))
  expect_true(all(summ[, "score"] %in% round(famscore_bdeu_byrow(tab, 1), 2)))
  expect_equal(summ[, "size"], rep(1, 4), ignore_attr = T)
})
