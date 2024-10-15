


files <- list.files(indir, ".rds")
indx <- grepl("barley", files)
files <- files[indx]


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
  MCMCchain <- res$MCMCchain

  N <- par$N
  r <- par$r

  # import bn
  set.seed(r)
  bn <- sim_load_bn(par)
  nlev <- vapply(bn, function(x) dim(x$prob)[1], integer(1))
  n    <- length(bn)

  # ground truth
  dag <- bnlearn::amat(bn)
  dmat <- bida:::descendants(dag)
  set.seed(r)
  pdo <- bida:::interv_probs_from_bn(bn, "bn")  # ground truth
  truetau <- matrix(vapply(pdo, bida:::JSD, numeric(1)), n, n)
  dindx <- diag(n) == 1
  tic[["ground truth"]] <- Sys.time()

  # draw data
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, N)
  tic[["simulate data"]] <- Sys.time()

  # compute support over unique dags
  burnin <- 1:200
  dags <- lapply(MCMCchain$traceadd$incidence[-burnin], as.matrix)
  tmp <- unique(dags)
  support <- bida:::rowsum_fast(rep(1/length(dags), length(dags)), dags, tmp)
  dags <- tmp

  # estimate intervention distributions ----
  ## compute support over parent sets
  ps <- bida::parent_support_from_dags(dags, support)
  tic[["compute parent support"]] <- Sys.time()

  get_size <- function(x) {
    dims <- x$counts$dim
    c(parents = length(dims)-1,
      parts = ifelse(is.null(x$partition), prod(dims[-1]), length(x$partition)))
  }

  ## compute mse of point-estimates (mean) of intervention distribution
  set.seed(r)
  mse <- tau <- parents <- parts <-  matrix(list(), n, n)
  runtime_fit <- runtime_mean <- matrix(0, n, n)
  cat("Start computing estimates for ", f, "\n")
  for (x in seq_len(n)) {
    cat(" Compute estimates for cause node x", x, "\n")
    for (y in seq_len(n)[-x]) {
      type <- ifelse(par$local_struct == "none", "cat", par$local_struct)
      pa   <- which(dag[, x] == 1) # true parents

      pairs <- list()
      tmp <- Sys.time()

      pairs$unknown <-list(
        unknown = bida::bida_pair(type, data, x, y,
                                  sets = ps$sets[[x]],
                                  support = ps$support[[x]],
                                  hyperpar = c(list(nlev = nlev), par)),
        full = bida::bida_pair("cat", data, x, y,
                               sets = ps$sets[[x]],
                               support = ps$support[[x]],
                               hyperpar = c(list(nlev = nlev), par)),
        known = bida::bida_pair(type, data, x, y,
                                sets = matrix(pa, nrow = 1),
                                support = 1,
                                hyperpar = c(list(nlev = nlev), par))
      )

      pdo_hat     <- lapply(pairs, bida::posterior_mean)
      mse[[x, y]]  <- vapply(pdo_hat, function(p) mean( (p-pdo[[x, y]])**2 ), numeric(1))
      tau[[x, y]]  <- vapply(pdo_hat, bida:::JSD, numeric(1))

      # number of conditioning variables - x + parents
      tmp <- lapply(pairs,
                    function(pair) do.call(rbind, lapply(pair$params, get_size))*pair$support)
      tmp <- vapply(tmp, colMeans, numeric(2))

      parents[[x, y]] <- tmp[1, ]
      parts[[x, y]]   <- tmp[2, ]
    }
	gc()

  }
  tic <- list("compute backdoor estimates" = Sys.time())

  out$mse_pdo <- colMeans(do.call(rbind, mse[!dindx]))
  out$mse_tau <- colMeans((do.call(rbind, tau[!dindx]) - truetau[!dindx])**2)
  out$parents <- colMeans(do.call(rbind, parents[!dindx]))
  out$parts   <- colMeans(do.call(rbind, parts[!dindx]))

  # edge probs and ranking ---
  edgep <- Reduce("+", Map("*", dags, support))
  arp   <- Reduce("+", Map("*", lapply(dags, bida:::descendants), support))
  taumat <- do.call(rbind, tau[!dindx])

  # compute average precision-recall
  compute_avgppv <- function(x, y) {
    indx <- order(x+runif(length(x))/1000, decreasing = TRUE)
    tp <- cumsum(y[indx])
    pp <- seq_along(x)
    mean((tp/pp)[y[indx] == 1])
  }
  eval_edge <- function(edgep, dag) {
    rates <- rowsum(edgep, dag)/tabulate(dag+1, 2)
    c(n = sum(edgep),
      fpr = rates[1],
      tpr = rates[2],
      avgppv = compute_avgppv(edgep, dag))
  }
  out$edgep <- eval_edge(edgep[!dindx], dag[!dindx])
  out$arp   <- eval_edge(arp[!dindx], dmat[!dindx])
  out$postau  <- apply((taumat>0)*1, 2, eval_edge, dag = dmat[!dindx])

  out$rank <- c(arp = compute_avgppv(arp[!dindx], dmat[!dindx]),
                apply(taumat, 2, compute_avgppv, y = dmat[!dindx]))
  topmat <- truetau > quantile(truetau[!dindx & truetau > 0], .8)
  out$ranktop <- c(arp = compute_avgppv(arp[!dindx], topmat[!dindx]),
                   apply(taumat, 2, compute_avgppv, y = topmat[!dindx]))
  tic <- list("evaluate" = Sys.time())

  out$par <- par
  out$tic <- diff(do.call(c, tic))
  out
}
