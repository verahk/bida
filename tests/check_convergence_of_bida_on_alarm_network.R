

bn <- readRDS("inst/data/alarm.rds")
nlev <- sapply(bnlearn::rbn(bn, 1), nlevels)
dags <- list(bnlearn::amat(bn))
dindx <- diag(length(bn)) == 1

pdo <- interv_probs_from_bn(bn, "exact")
pair  <- rep.int(seq_along(dindx)[!dindx], lengths(pdo[!dindx]))
truth <- list(pdo = pdo[!dindx],
              jsd = vapply(pdo[!dindx], avg_jsd_array, numeric(1)))

samplesizes <- round(10**seq(2, 4, length.out = 5))
rmse <- list(pdo = list(), jsd = list())
N <- 100

rmse <- list(pdo = matrix(nrow = length(truth$pdo), ncol = length(samplesizes)),
             jsd = matrix(nrow = length(truth$pdo), ncol = length(samplesizes)))

#for (x in seq_len(n)) for(y in seq_len(n)[-x]) stopifnot(all(dim(tmp[[x, y]]) == dim(pdo[[x, y]])))
for (i in seq_along(samplesizes)) {

  N <- samplesizes[i]
  data <- sample_data_from_bn(bn, N)
  bp <- bida_sample(data, dags, type = "categorical", adjset = "o_min",
                    params = list(ess = 1, nlev = nlev, maxconf = Inf))

  pdo <- posterior_mean(bp)
  rmse$pdo[, i]  <- mapply(function(x, y) sqrt(mean((x-y)**2)),
                           x = pdo[!dindx],
                           y = truth$pdo)

  jsd <- posterior_mean(bp, n = 10**3, ace_funs = list(avg_jsd_array))
  rmse$jsd[, i] <- (unlist(jsd[!dindx])-truth$jsd)**2
}

matplot(t(rmse$pdo), type = "l", col = "grey",
        xlab = "sample size", ylab = "RMSE", main = "Posterior mean, intervention probs")
matplot(t(rmse$jsd), type = "l", col = "grey",
        xlab = "sample size", ylab = "RMSE", main = "Posterior mean, jsd"))
