

r <- 2
alpha <- c(1/seq(10, 1), 2, 3)
plot(alpha, lgamma(r*alpha), type = "l")
lines(alpha, r*lgamma(alpha), col = "red")

levels <- list(0:1, 0:1, 0:1)
conf   <- expand.grid(levels)
stride <- c(1, cumprod(lengths(levels[-length(levels)])))
r <- 2
q <- nrow(conf)
m <- cbind(c(rep(0, 4), rep(10, 4)), c(rep(0, 4), c(10, 10, 100, 100)))

# no partitioning
tab <- m
scores <- famscore_bdeu_byrow(tab, 1, r, q, s = 1)
cbind(tab, scores, sum(scores))

# fully collapsed
group <- rep(0, each = q)
tab <- rowsum(m, group)
scores <- famscore_bdeu_byrow(tab, 1, r, q, s = tabulate(group+1))
cbind(tab, scores, sum(scores))


# split on one variable
for (i in seq_along(levels)) {
  group <- conf[, i]
  tab <- rowsum(m, group)
  scores <- famscore_bdeu_byrow(tab, 1, r, q, s = tabulate(group+1))
  print(cbind(tab, scores, sum(scores)))
}


# collapse two rows
for (i in seq_along(levels)) {
  group <- as.matrix(conf[, -i])%*%stride[-i]
  tab <- rowsum(m, group)
  scores <- famscore_bdeu_byrow(tab, 1, r, q, s = tabulate(match(group, unique(group))))
  print(cbind(tab, scores, sum(scores)))
}

