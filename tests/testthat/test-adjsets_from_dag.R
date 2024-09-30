

test_that("adjset bida and pcalg return same sets", {

  ## Check that adjustment sets returned by bida:::adjset_from_dag
  ## coincide with those from pcalg::optAdjSet (o-set) and pcalg::adjustment (minimal sets)

  set.seed(007)
  adjsets <- c("o", "o_min", "pa_min", "pa", "anc")
  n <- 10
  ngraphs <- 1
  verbose <- F

  ## o-set
  for (g in 1:ngraphs){

    dag <- as(pcalg::randDAG(n, 3, weighted = FALSE), "matrix")
    colnames(dag) <- rownames(dag) <- paste0("X", seq_len(n))
    dmat  <- descendants(dag)
    tdag <- t(dag)

    # identify all adjustment set
    sets <- adjsets_from_dag(adjsets, dag, dmat)

    #pcalg::plot(as(dag, "graphNEL"))
    for (x in 1:n){
      for (y in seq_len(n)[-x]){
        if (dmat[x, y] == 0) {
          exp <- rep(list(y), 3)
          expect_equal(sets[x, y, c("pa_min", "o", "o_min")], exp, ignore_attr = T)
          expect_equal(sets[[x, y, "pa"]], which(dag[, x] == 1), ignore_attr = T)
        } else {
          if (verbose) cat("\ng = ", g, "x = ", x, "y = ", y)

          # o-set
          expect_equal(sort(sets[[x, y, "o"]]),
                       sort(pcalg::optAdjSet(tdag, x, y)))

          # minimal sets
          o  <- sort(sets[[x, y, "o_min"]])
          pa <- sort(sets[[x, y, "pa_min"]])

          min <- pcalg::adjustment(tdag, "dag", x, y, "minimal")

          # loop through all minimal sets and compare with o/pa
          o_equal <- pa_equal <- FALSE
          for (z in min) {
            z <- sort(z)
            if (!o_equal)  o_equal <- length(o) == length(z) && all(o == z)
            if (!pa_equal) pa_equal <- length(pa) == length(z) && all(pa == z)
            if (o_equal && pa_equal) break
          }
          expect_equal(o_equal, TRUE)
          expect_equal(pa_equal, TRUE)
        }
      }
    }
  }
})


test_that("function replace large sets with smaller", {

  dag <- rbind(U   = c(0, 0, 0, 0, 1),
               Z1  = c(0, 0, 0, 1, 1),
               Z2  = c(0 ,0, 0, 1, 1),
               X   = c(0 ,0, 0, 0, 1),
               Y   = c(0, 0, 0, 0, 0))
  colnames(dag) <- rownames(dag)
  #Rgraphviz::plot(as(dag, "graphNEL"))

  x <- 4
  y <- 5
  adjsets <- c("o", "o_min", "pa_min")
  nlev <- rep(3, ncol(dag))

  # no checksize
  sets <- adjsets_from_dag(adjsets, dag, xvars = x, yvars = y)
  expect_equal(sets$pa_min, 2:3, ignore_attr = T)
  expect_equal(sets$o_min, 2:3, ignore_attr = T)
  expect_equal(sets$o, 1:3, ignore_attr = T)

  # maximum conf of adjustment set is 9 --> replace o with o_min
  maxconf <- 3**4/outer(nlev, nlev)
  sets <- adjsets_from_dag(adjsets, dag, xvars = x, yvars = y,
                           checksize = function(x, y, z) length(z) <= 1 || prod(nlev[z]) <= maxconf[x, y])
  expect_equal(sets$pa_min, 2:3, ignore_attr = T)
  expect_equal(sets$o_min, 2:3, ignore_attr = T)
  expect_equal(sets$o, 2:3, ignore_attr = T)

  # maximum conf of adjustment set is 3 --> replace all sets with NA
  maxconf <- 3**3/outer(nlev, nlev)
  sets <- adjsets_from_dag(adjsets, dag, xvars = x, yvars = y,
                           checksize = function(x, y, z) length(z) <= 1 || prod(nlev[z]) <= maxconf[x, y])
  expect_true(is.na(sets$pa_min))
  expect_true(is.na(sets$o_min))
  expect_true(is.na(sets$o))
})


