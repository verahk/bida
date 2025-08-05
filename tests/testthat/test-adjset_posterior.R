
test_that("function returns expected output in example dag", {

  dag <- rbind(U   = c(0, 0, 0, 0, 1),
               Z1  = c(0, 0, 0, 1, 1),
               Z2  = c(0 ,0, 0, 1, 1),
               X   = c(0 ,0, 0, 0, 1),
               Y   = c(0, 0, 0, 0, 0))
  colnames(dag) <- rownames(dag)
  #Rgraphviz::plot(as(dag, "graphNEL"))
  dag2 <- dag
  dag2["X", ] <- 0

  pdags <- list(dags = list(dag, dag2),
                p = c(.5, .5))

  x <- 4
  y <- 5

  # minimal parents
  adjset <- "pa_min"
  expected <- list(sets = list(y, c(2, 3)), p = c(.5, .5), zeroprob = .5)
  attr(expected, "adjset") <- adjset
  expect_equal(adjset_support(pdags, x, y, adjset), expected)

  # o-set
  adjset <- "o"
  expected <- list(sets = list(y, c(1, 2, 3)), p = c(.5, .5), zeroprob = .5)
  attr(expected, "adjset") <- adjset
  expect_equal(adjset_support(pdags, x, y, adjset), expected)

  # minimal o-set
  adjset <- "o_min"
  expected <- list(sets = list(y, c(2, 3)), p = c(.5, .5), zeroprob = .5)
  attr(expected, "adjset") <- adjset
  expect_equal(adjset_support(pdags, x, y, adjset), expected)

  # parent set
  adjset <- "pa"
  expected <- list(sets = list(c(2, 3)), p = 1)
  attr(expected, "adjset") <- adjset
  expect_equal(adjset_support(pdags, x, y, adjset), expected)

})


