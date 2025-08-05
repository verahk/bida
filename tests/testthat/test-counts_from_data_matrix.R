
test_that("function gives same res as table", {

  # sample a set of categorical levels with some levels unobserved
  data <- replicate(10, sample(0:2, size = 2, TRUE))
  df <- data.frame(apply(data, 2, factor, levels = 0:2, simplify = FALSE))

  expect_equal(counts_from_data_matrix(data, rep(3, 10)),
               unclass(unname(table(df))))
})

