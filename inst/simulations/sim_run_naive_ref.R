par <- list(local_structure = NA,
            init = NA,
            sample = NA,
            ess = 1,
            edgepf = NA,
            hardlimit = NA,
            r = 1:30)

# add params controlling random-cpt generation
par$n <- c(20)
par$maxdepth <- c(1)
par$k <- c(4)

pargrid <- expand.grid(par, stringsAsFactors = FALSE)

params_to_filename <- function(par) {
  tmp <- sprintf("n%s_k%s_depth%s_%s_%s_%s_ess%s_epf%s_N%s_r%02.0f.rds",
                 par$n, par$k, par$maxdepth*100, par$init, par$local_struct, par$sample, par$ess, par$edgepf, par$N, par$r)
  stopifnot(length(tmp) == 1)  # fails if any argument is NULL
  tmp
}

for (r in 1:nrow(pargrid)) {
  print(unlist(pargrid[r, ]))
  sim_run(pargrid[r, ], bn, outdir)
 # profvis::profvis(sim_run(pargrid[r, ], bn, outdir))
}



sim_run <- function(par, bn, outdir, verbose = FALSE) {


  out <- list()
  r <- par$r

  cat("Load bn\n")
  set.seed(r)
  bn <- sim_load_bn(par)
  nlev <- vapply(bn, function(x) dim(x$prob)[1], integer(1))
  n    <- length(bn)

  cat("Compute ground truth\n")
  dag <- bnlearn::amat(bn)
  dmat <- bida:::descendants(dag)
  set.seed(r)
  pdo <- bida:::interv_probs_from_bn(bn, "bn")  # ground truth
  truetau <- matrix(vapply(pdo, bida:::JSD, numeric(1)), n, n)
  dindx <- diag(n) == 1

  # compute average precision-recall
  compute_avgppv <- function(x, y) {
    indx <- order(x+runif(length(x))/1000, decreasing = TRUE)
    tp <- cumsum(y[indx])
    pp <- seq_along(x)
    mean((tp/pp)[y[indx] == 1])
  }

  # draw data
  for (N in c(300, 1000, 3000)) {
    set.seed(N+r)
    data <- bida:::sample_data_from_bn(bn, N)


    ## marginal probs -----
    out <- list()
    par$N <- N
    cat("Estimate intervention distrib for sample size", N, "\n")

    for (ref in c("marg", "cond")) {
      par$local_struct <- ref

      # compute estimates and mse of interv probs
      mse <- tau <- parents <- parts <- matrix(list(), n, n)
      if (ref == "marg") {
        for (y in seq_len(n)) {
          bdeu<- bida:::bida_bdeu(data, y, integer(0), ess = 1, nlev)
          pdo_hat   <- c(bida:::backdoor_mean(bdeu, 1))
          mse[-y, y] <- lapply(pdo[-y, y], function(p) mean((p-pdo_hat)**2))
          tau[-y, y] <- list(0)
        }
      } else {
        for (x in seq_len(n)) {
          for (y in seq_len(n)[-n]) {
            bdeu <- bida:::bida_bdeu(data, y, x, ess = 1, nlev)
            pdo_hat   <- bida:::backdoor_mean(bdeu)
            mse[[x, y]] <- mean((pdo_hat-pdo[[x, y]])**2)
            tau[[x, y]] <- bida:::JSD(pdo_hat)
          }
        }
      }

      # evaluate
      out$mse_pdo <- c(unknown = mean(unlist(mse[!dindx])))
      tmp <- compute_avgppv(unlist(tau[!dindx]), dmat[!dindx])
      out$rank    <-c(arp = tmp, unknown = tmp)
      topmat <- truetau > quantile(truetau[!dindx & truetau > 0], .8)
      tmp <- compute_avgppv(unlist(tau[!dindx]), topmat[!dindx])
      out$ranktop <- c(arp = tmp, unknown = tmp)

      # store to file
      out$par <- par
      filename <- params_to_filename(par)
      saveRDS(out, paste0(outdir, filename))
    }
  }
}

stop()
