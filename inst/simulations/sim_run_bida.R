


files <- list.files(outdir, ".rds")
indx <- grepl("n20_k4.*depth50.*ptree.*logN.*", files)
print(sum(indx))
#print(files[indx])


files <- list.files(indir, ".rds")
files <- files[grepl("epf2", files) & grepl("water2|insurance|asia", files)]
files <- files[! (grepl("N100_", files) | grepl("r2", files)) ]
files <- files[!file.exists(paste0(outdir, files))]
cat("Number of files:", length(files), "\n")
pargrid <- data.frame(file = files,
                      indir = indir)

params_to_filename <- function(par) par$file

sim_run <- function(par, verbose = FALSE) {

  indir <- par$indir
  f     <- par$file

  out <- list()
  tic <- list(start = Sys.time())

  # read MCMC-results from file
  filepath <- paste0(indir, f)
  res <- readRDS(filepath)
  par <- res$par
  lookup <- res$lookup

  N <- par$N
  r <- par$r

  # import bn
  set.seed(r)
  if (verbose) cat("Load bn\n")
  bn <- sim_load_bn(par)
  nlev <- vapply(bn, function(x) dim(x$prob)[1], integer(1))
  n    <- length(bn)
  tic[["load bn"]] <- Sys.time()

  if (verbose) cat("Compute ground truth bn\n")
  dag <- bnlearn::amat(bn)
  dmat <- bida:::descendants(dag)
  mode(dag) <- mode(dmat) <- "integer"
  set.seed(r)
  pdo <- bida:::interv_probs_from_bn(bn, "bn")  # ground truth
  tau <- matrix(vapply(pdo, bida:::JSD, numeric(1)), n, n)
  dindx <- diag(n) == 1
  tic[["ground truth"]] <- Sys.time()

  if (verbose) cat("Simulate data\n")
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, N)
  tic[["simulate data"]] <- Sys.time()

  ## compute support over parent sets ---
  if (verbose) cat("Compute parent support\n")
  MCMCchain <- res$MCMCchain
  burnin <- 1:200
  dags <- lapply(MCMCchain$traceadd$incidence[-burnin], as.matrix)
  tmp <- unique(dags)
  support <- bida:::rowsum_fast(rep(1/length(dags), length(dags)), dags, tmp)
  dags <- tmp
  ps <- bida::parent_support_from_dags(dags, support)
  out$parents <- mapply(function(m, w) sum(rowSums(m>0, na.rm = T)*w), ps$sets, ps$support)
  tic[["compute parent support"]] <- Sys.time()

  # estimate intervention distributions ----
  if (verbose) cat("Estimate intervention distributions\n")
  # helper functions
  compute_mse <- function(x, y) mean( (x-y)**2 )
  compute_avg_number_of_params <- function(pair){
    # compute avg number of params in cpt P(Y|X, Z) for all adjsets Z
    tmp <- vapply(pair$params, function(counts) prod(dim(counts)), double(1))
    sum(tmp*pair$support)
  }
  mse <- tau_hat <- params <- runtime <- list()
  for (name in c("unknown", "full", "known")) {
    if (name == "known") {
      ps <- bida:::parent_support_from_dags(list(dag), 1)
    }

    tictic <- Sys.time()
    type <- ifelse(par$local_struct == "none" || name == "full", "cat", par$local_struct)
    bida    <- bida:::bida(ps, data, type, list(nlev = nlev, ess = par$ess), lookup = NULL)
    pdo_hat <- bida:::posterior_mean(bida)
    runtime[[name]] <- Sys.time()-tictic

    mse[[name]]     <- mapply(compute_mse, pdo_hat[!dindx], pdo[!dindx])
    tau_hat[[name]] <- vapply(pdo_hat[!dindx], bida:::JSD, numeric(1))

    tmp <- vapply(bida[!dindx], compute_avg_number_of_params, numeric(1))
    params[[name]] <- mean(tmp)
  }

  tic[["compute interv probs"]] <- Sys.time()

  out$mse_pdo <- colMeans(do.call(cbind, mse))
  out$tau    <- cbind(do.call(cbind, tau_hat))
  out$mse_tau <- colMeans((out$tau - tau[!dindx])**2)
  out$params <-  unlist(params)
  out$runtime <- unlist(runtime)

  # edge probs and ranking ---
  edgep <- Reduce("+", Map("*", dags, support))
  arp   <- Reduce("+", Map("*", lapply(dags, bida:::descendants), support))

  # compute average precision-recall
  compute_avgppv <- function(x, y) {
    indx <- order(x+runif(length(x))/1000, decreasing = TRUE)
    tp <- cumsum(y[indx])
    pp <- seq_along(x)
    mean((tp/pp)[y[indx] == 1])
  }
  eval_edge <- function(edgep, dag) {
    rates <- bida:::rowsum_fast(edgep, dag, 0:1)/tabulate(dag+1, 2)
    c(n = sum(edgep)/sum(dag),
      fpr = rates[1],
      tpr = rates[2],
      avgppv = compute_avgppv(edgep, dag))
  }
  out$edgep <- eval_edge(edgep[!dindx], dag[!dindx])
  out$arp   <- eval_edge(arp[!dindx], dmat[!dindx])
  out$rank <- c(arp = compute_avgppv(arp[!dindx], dmat[!dindx]),
                apply(out$tau, 2, compute_avgppv, y = dmat[!dindx]))
  topmat <- tau > quantile(tau[!dindx & tau > 0], .8)
  out$ranktop <- c(arp = compute_avgppv(arp[!dindx], topmat[!dindx]),
                   apply(out$tau, 2, compute_avgppv, y = topmat[!dindx]))
  tic[["eval"]] <- Sys.time()

  out$par <- par
  out$tic <- diff(do.call(c, tic))
  out
}
