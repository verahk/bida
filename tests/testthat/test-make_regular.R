test_that("multiplication works", {
  P <- list(0:3)
  nlev <- c(2, 2)
  obj <- make_regular(P, nlev)
  exp <- as.list(0:3) 
  expect_true(all(match(obj, exp, 0) > 0))
})
