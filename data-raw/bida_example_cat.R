## code to prepare `bida_example` dataset goes here

# specify DAG row-wise:
dag <- rbind(Z1  = c(0, 0, 0, 1, 0, 0, 0),
             Z2  = c(0, 0, 0, 1, 0, 0, 0),
             L   = c(0, 1, 0, 0, 1, 0, 0),
             X   = c(0, 0, 0, 0, 1, 0, 0),
             M   = c(0, 0, 0, 0, 0, 1, 0),
             Y   = c(0, 0, 0, 0, 0, 0, 0),
             U   = c(0, 0, 0, 0, 0, 1, 0))
colnames(dag) <- rownames(dag)

graphNEL <- as(dag, "graphNEL")

# define bn
set.seed(42)
n <- ncol(dag)
nlev <- rep(3, n)
cpts <- rand_dist(dag, "cat", list(nlev = nlev, alpha = 1))
bn   <- custom_bn(dag, cpts)

# sample data
df <- bnlearn::rbn(bn, 300)
data <- vapply(df, as.integer, integer(nrow(df)))-1L

# run partition MCMC
scorepar  <- BiDAG::scoreparameters("bdecat", data = data.frame(lapply(df, factor, exclude = NULL)))
MCMCchain <- BiDAG:::sampleBN(scorepar, "partition")

# collect objects in list and write to file
bida_example_cat <- list(dag = dag,
                         nlev = nlev,
                         cpts = cpts,
                         bn = bn,
                         data = data,
                         partitionMCMC = MCMCchain)
usethis::use_data(bida_example_cat, overwrite = TRUE)
