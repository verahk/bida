

r <- 3
q <- 11
x <- c(.5, 1:10)

plot(x, 2*lgamma(x))
lines(x, lgamma(2*x))

counts <- matrix(x, q, r)

score1 <- bida:::famscore_bdeu_byrow(counts, 1, r, q, s = 1)
score2 <- bida:::famscore_bdeu_byrow(counts, 1, r, q, s = 2)
plot(x, score1)
lines(x, score2, col = "red")

# two rows with zero counts
s <- 2
m <- matrix(0, s, r)
score_rows <- bida:::famscore_bdeu_byrow(m, 1, r, q, s = 1)
bida:::famscore_bdeu_byrow(matrix(0, 1, r), 1, r, q, s = s)
score_collapsed <- bida:::famscore_bdeu_byrow(matrix(colSums(m), 1), 1, r, q, s = s)
