

test_that("score equivalence of bdeu score", {

  dag1 <- rbind(z = c(0, 1, 1),
                x = c(0, 0, 1),
                y = c(0, 0, 0))
  dag2 <- rbind(z = c(0, 0, 1),
                x = c(1, 0, 1),
                y = c(0, 0, 0))
  dag3 <- rbind(z = c(0, 0, 0),
                x = c(1, 0, 0),
                y = c(1, 1, 0))

  params <- list(ess = 1,
                 nlev = c(2, 2, 2))

  # test score equivalence standard DAG
  data <- sapply(params$nlev, sample, size = 100, replace = T)-1
  expect_equal(score_dag(data, dag1,  type = "cat", params),
               score_dag(data, dag2, type = "cat", params))
  expect_equal(score_dag(data, dag1, type = "cat", params),
               score_dag(data, dag3, type = "cat", params))

})
