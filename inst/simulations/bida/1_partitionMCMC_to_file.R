# WHAT: run partition MCMC
# WHY: intermediate step in bida-simulations, with long runtime.
# HOW:
# For each network, samplesize and iteration do:
# 1. read in netowork from file
# 2. run partitionMCMC
# 3. save partitionMCMC-object to file


library(foreach)
library(bida)

getwd()
par <- expand.grid(bnname = c("child",
                              "alarm",
                              "win95pts"),
                   N = c(300, 1000, 3000),
                   r = 1:30)

# set up cluster
ncores <- 4
cl <- parallel::makeCluster(ncores, type="SOCK", outfile = "")
doSNOW::registerDoSNOW(cl)
on.exit(parallel::stopCluster(cl))

# define function
partitionMCMC_to_file <- function(bn, N, r) {

  bn <- readRDS(paste0("../bida/inst/data/", bnname, ".rds"))
  tag <-sprintf("%s_N%0.0f_r%02.0f", bnname, N, r)
  filepath <- paste0("./", tag, ".rds")
  if (file.exists(filepath)) return(NULL)
  cat("\n", format(Sys.time(), "%a %b %d %H:%M:%S"), "tag:", tag)

  set.seed(r+N)
  tic <- Sys.time()
  data <- bida:::sample_data_from_bn(bn, N)              # sample data, until at least 2 levels observed
  colnames(data) <- names(bn)
  tic <- c(tic, data = Sys.time())

  # run MCMC
  df <- data.frame(lapply(asplit(data, 2), factor, exclude = NULL))   # remove unused levels
  my_score <- BiDAG::scoreparameters(scoretype = "bdecat",
                                     data = df,
                                     bdecatpar = list(chi=1,edgepf=1))

  partfit <- BiDAG::partitionMCMC(scorepar = my_score,
                                  startspace = NULL)
  tic <- c(tic, partitionMCMC = Sys.time())
  attr(partfit, "toc") <- diff(tic)
  attr(partfit, "Sys.info") <- Sys.info()

  saveRDS(partfit, filepath)
}

# test
if (FALSE) {
  bnname <- "asia"
  N <- 300
  r <- 0

  wd <- getwd()
  setwd("../partitionMCMC")
  partitionMCMC_to_file(bnname, N, r)

  filepath <- sprintf("%s_N%0.0f_r%02.0f.rds", bnname, N, r)
  test <- readRDS(filepath)
  names(test)
  attr(test, "toc")
  test$info
  file.remove(filepath)
  setwd(wd)
}

# run
if (ncores > 0) {
  foreach (bnname = par$bnname,
           N = par$N,
           r = par$r) %dopar% partitionMCMC_to_file(bnname, N, r)
} else {
  for (r in 1:nrow(par)) partitionMCMC_to_file(par$bnname[r], par$N[r], par$r[r])
}


