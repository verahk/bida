
# What: Change distribution of bayesian network water.rds
# Why:  Original bn has source nodes with only one outcome.
#       Need data with multiple levels to apply BiDAG and pcalg structure learning functions
# How:  Generate a copy of bayesian network water.rds with uniform distribution for every source node

bn  <- readRDS("./inst/data/water.rds")
dag <- bnlearn::amat(bn)
dist  <- lapply(bn, "[[", "prob")

# all nodes with 0-1 probs are source nodes
indx <- sapply(dist, function(x) all(x %in% 0:1))
par  <- lapply(bn[indx], "[[", "parents")
all(lengths(par) == 0)

# generate copy of distribution with uniform distrib for source nodes
new_dist <- dist
for (i in which(indx)) {
  new_dist[[i]][] <- 1/length(dist[[i]])
}
new_dist[indx]

# generate new bn object
g <- bnlearn::empty.graph(names(bn))
bnlearn::amat(g) <- dag
new_bn <- bnlearn::custom.fit(g, new_dist)


# test
df <- bnlearn::rbn(new_bn, 10**4)
lapply(df[indx], table)

saveRDS(new_bn, file = "./inst/data/water2.rds")
