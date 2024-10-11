
test_that("score function is score equivalent", {
  data <- data.frame(replicate(3, sample(1:3, 10^3, replace = TRUE)))
  names(data) <- paste0("X", seq_along(data))

  # X1 -> X2 -> X3
  m1 <- matrix(table(data[, 1]),   ncol = 3, byrow = TRUE) # X1
  m2 <- matrix(table(data[, 2:1]), ncol = 3, byrow = TRUE) # X2 | X1,
  m3 <- matrix(table(data[, 3:2]), ncol = 3, byrow = TRUE) # X3 | X2
  sc1 <- famscore_bdeu(m1) + famscore_bdeu(m2) + famscore_bdeu(m3)

  # X1 <- X2 <- X3
  m1 <- matrix(table(data[, 1:2]),   ncol = 3, byrow = TRUE)
  m2 <- matrix(table(data[, 2:3]), ncol = 3, byrow = TRUE)
  m3 <- matrix(table(data[, 3]), ncol = 3, byrow = TRUE)
  sc2 <- famscore_bdeu(m1) + famscore_bdeu(m2) + famscore_bdeu(m3)

  # X1 <- X2 -> X3
  m1 <- matrix(table(data[, 1:2]),   ncol = 3, byrow = TRUE)
  m2 <- matrix(table(data[, 3:2]), ncol = 3, byrow = TRUE)
  m3 <- matrix(table(data[, 2]), ncol = 3, byrow = TRUE)
  sc3 <- famscore_bdeu(m1) + famscore_bdeu(m2) + famscore_bdeu(m3)

  expect_equal(sc1, sc2)
  expect_equal(sc1, sc3)
})

test_that("famscore_bdeu_byrow correctly scorse sparse matrix", {

  score_sparse_matrix <- function(n, ess, r, q, s) {
    m <- matrix(0, length(s), r)
    m[1] <- n
    famscore_bdeu_byrow(m, ess, r, q, s)
  }
  r <- 4
  q <- r**3
  full <- score_sparse_matrix(1, 1, r, q, rep(1, q))
  part <- score_sparse_matrix(1, 1, r, q, rep(q/r**2, r))
  collapse <- score_sparse_matrix(1, 1, r, q, rep(q, 1))

  expect_equal(sum(full), sum(part))
  expect_equal(sum(full), sum(part))
  expect_equal(sum(full), sum(collapse))
  expect_equal(sum(full), famscore_bdeu_1row(c(1, rep(0, r-1)), 1, r))
})
