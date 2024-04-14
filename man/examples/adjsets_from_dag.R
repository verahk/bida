
# specify DAG row-wise:
dag <- rbind(Z1  = c(0, 0, 0, 1, 0, 0, 0),
             Z2  = c(0, 0, 0, 1, 0, 0, 0),
             L   = c(0, 1, 0, 0, 1, 0, 0),
             X   = c(0, 0, 0, 0, 1, 0, 0),
             M   = c(0, 0, 0, 0, 0, 1, 0),
             Y   = c(0, 0, 0, 0, 0, 0, 0),
             U   = c(0, 0, 0, 0, 0, 1, 0))
colnames(dag) <- rownames(dag)
g <- graph::graphAM(dag, edgemode = "directed")


# compute adjustment set w.r.t. X and Y
x <- 4
y <- 6
adjsets <- c("max", "pa","pa_min", "o", "o_min")
sets <- bida:::adjsets_from_dag(adjsets, dag, xvars = x, yvars = y)

# plot each adjustment set
color <- c("grey", "lightblue",  "blue", "pink", "red")
names(color) <- adjsets
par(mfrow = c(1, 5),
    mar = c(5, 1, .1, .1))

for (a in names(sets)) {
  z <- sets[[a]]
  tmp <- rep(color[a], length(z))
  names(tmp) <- colnames(dag)[z]
  Rgraphviz::plot(g,
                  nodeAttrs = list(fillcolor = tmp),
                  main = a)
}



# define function for testing adjustment set size
test_adjset_size_cat <- function(adjset, x, y, z, testarg) {
  is.na(z) || prod(testarg$nlev[c(x, y, z)]) < testarg$maxconf
}

adjsets <- c("pa", "o", "o_min", "pa_min")
adjsets_from_dag(adjsets, dag)


checksize <- list()
checksize$fun <- check_adjset_size_cat <- function(adjset, x, y, z, checksize) {
  is.na(z) || prod(checksize$nlev[c(x, y, z)]) < checksize$maxconf
}
checksize$maxconf <- 10
