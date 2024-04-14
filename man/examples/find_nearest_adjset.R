# specify DAG by children/row-wise:
dag <- rbind(L = c(0, 1, 0, 1, 0),
             Z = c(0, 0, 1, 0, 0),
             X = c(0, 0, 0, 1, 1),
             Y = c(0, 0, 0, 0, 1),
             C = c(0, 0, 0, 0, 0))
colnames(dag) <- rownames(dag)
g <- as(dag, "graphNEL")
Rgraphviz::plot(g)


n <- ncol(dag)
dmat <- bida::descendants(dag)

X <- 3
Y <- 4

# nodes reachable from X via all active paths
reachable <- find_nodes_reachable_via_ancestors(dag, dmat, X-1, A = rep(TRUE, n), Z = 0L)
colnames(dag)[reachable]


# nodes reachable from X via active ancestor paths
A <- dmat[, X] == 1 | dmat[, Y] == 1
reachable <- find_nodes_reachable_via_ancestors(dag, dmat, X-1, A, Z = 0L)
colnames(dag)[reachable]

# find adjustment sets w.r.t. X and Y
G <- dag
G[X, ] <- 0                         # backdoor graph
Z0 <- which(A & dmat[X, ] == 0)     # ancestor excluding descendants
pa <- find_nearest_adjset(G, dmat, X-1, A, Z0-1) +1
o  <-  find_nearest_adjset(G, dmat, Y-1, A, Z0-1) +1

