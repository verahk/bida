
# WHAT: check convergence of parameters estimated with functions bida:::backdoor_params_cat
# WHY: use random input to check functions works as expected
# HOW:
# - draw a random BN (ground truth)
# - repeatedly sample data from the BN and estimate params
# - plot mean squeared errors against increasing sample size

library(Rgraphviz) # for plotting
library(ggplot2)
library(dplyr)
library(bida)

# DRAW A RANDOM DAG ----
# - specify a BN over three variables
# - compute intervention probs and causal effects
#set.seed(007)
dag <- rbind(c(0, 1, 0, 0),
             c(0, 0, 0, 0),
             c(1, 1, 0, 0),
             c(1, 1, 0, 0))
colnames(dag) <- rownames(dag) <- c("X", "Y", "Z1", "Z2")
dmat <- bida:::descendants(dag)
plot(as(dag, "graphNEL"))


# draw a distribution over the DAG
n <- ncol(dag)
nlev <- c(2, 3, 4, 5)
cpts <- rand_cpt_arrays(dag, nlev, alpha= 1)
bn   <- bn_from_cpt_arrays(cpts)
data <- sample_data_from_bn(bn, 10**3)

py   <- lapply(seq_len(n), function(y) tabulate((data+1)[, y], nlev[y])/10**3)

# compute true causal params
pdo <- interv_probs_from_cpt_arrays(cpts, "exact")
tau <- list(jsd = array(vapply(pdo, avg_jsd_array, numeric(1)), dim(pdo)),
            abs = array(vapply(pdo, avg_abs_array, numeric(1)), dim(pdo)),
            dev = array(Map(
              function(x, y) {if (dmat[x, y]==0) return(0); abs_dev(pdo[[x, y]], M = py[[y]])},
                              x = rep(seq_len(n), times = n),
                              y = rep(seq_len(n), each = n)),
                        dim(pdo)))




# ESTIMATE PARAMS ----
# - repeatedly sample data from the BN, estimate backdoor params, and compute mean squared errors.

pairs <- list(adjust = list(x = 1, y = 2,  z = 3:4),
              cond = list(x = 3, y = 2, z = integer(0)),
              marg = list(x = 3, y = 4, z = 4))
simpar <- list("class" = c("array", "sparse"),
               "fun" = c("posterior_mean", "posterior_sample", "jsd", "abs", "dev_marg"),
               "samplesize" = round(10**seq(1, 4, length.out = 10)),
               "iter" = 1:30)



for (p in names(pairs)) {
  x <- pairs[[p]]$x
  y <- pairs[[p]]$y
  z <- pairs[[p]]$z
  mse <- array(dim = lengths(simpar), dimnames = simpar)

  cat("\npair:", p)
  print(pairs[[p]])

  for (NN in seq_along(simpar$samplesize)) {
    N <- simpar$samplesize[NN]
    for (r in simpar$iter){

      # sample data and compute params
      #set.seed(r+N)
      data <- sample_data_from_bn(bn, N)#bida:::sample_data_from_cpts(cpts, N)
      ay   <- tabulate(data[, y] + 1, nlev[y]) + 1/nlev[y]
      fit <- list(array  = bida:::backdoor_params_cat(data, x, y, z, nlev, min_sparse = Inf),
                  sparse = bida:::backdoor_params_cat(data, x, y, z, nlev, min_sparse = 0))

      for (cl in seq_along(fit)) {

        # analytical mean
        tmp <- posterior_mean.backdoor_params_cat(fit[[cl]], ess = 1, kx = nlev[x])
        mse[cl, 1, NN, r] <- mean((tmp-pdo[[x, y]])**2)

        smpl <- posterior_sample.backdoor_params_cat(fit[[cl]], samplesize = 10**3, ess = 1, kx = nlev[x])

        # mean of sample
        tmp <- rowMeans(smpl, dims = 2)
        mse[cl, 2, NN, r] <- mean((tmp-pdo[[x, y]])**2)

        # causal effects
        tmp <- mean(bida:::avg_jsd_array(smpl))
        mse[cl, 3, NN, r] <- mean((tmp-tau$jsd[[x, y]])**2)

        tmp <- mean(bida:::avg_abs_array(smpl))
        mse[cl, 4, NN, r] <- mean((tmp-tau$abs[[x, y]])**2)


        py  <- t(bida:::rDirichlet(10**3, ay, nlev[y]))
        tmp <- mean(abs_dev(smpl, k = dim(smpl), M = py))
        mse[cl, 5, NN, r] <- mean((tmp-tau$dev[[x, y]])**2)
      }
    }
  }

  # PLOT MSE ----
  df <- reshape2::melt(sqrt(mse)) %>%
    group_by(samplesize, fun, class) %>%
    summarize(min = min(value),
              max = max(value),
              mean = mean(value), .groups = "drop")

  df %>%
    select(samplesize, class, fun, mean) %>%
    tidyr::pivot_wider(names_from = "samplesize", values_from = "mean") %>%
    arrange(class, fun)

  title <- sprintf("%s : %s->%s (adjustment set: %s)",
                   p, colnames(dag)[x], colnames(dag)[y], paste(colnames(dag)[z], collapse = ","))
  plot <- ggplot(df, aes(samplesize, ymin = min, ymax = max, y = mean, color = class, fill = class)) +
    facet_grid(fun~class, scales = "free") +
    scale_x_log10() +
    #scale_y_log10() +
    geom_ribbon(alpha = .25) +
    geom_line() +
    geom_point() +
    ylab("rmse") +
    ggtitle(title)

  print(plot)
}




# PROFILING ----
x <- pairs[[1]]$x
y <- pairs[[1]]$y
z <- pairs[[1]]$z
fit <- list(array  = bida:::backdoor_params_cat(data, x, y, z, nlev, min_sparse = Inf),
            sparse = bida:::backdoor_params_cat(data, x, y, z, nlev, min_sparse = 0))


pryr::object_size(fit$array)
pryr::object_size(fit$sparse)

microbenchmark::microbenchmark(
  bida:::backdoor_params_cat(data, x, y, z, nlev, min_sparse = Inf),
  bida:::backdoor_params_cat(data, x, y, z, nlev, min_sparse = 0)
)


microbenchmark::microbenchmark(
  posterior_sample.backdoor_params_cat(fit$array, 10**3, ess = 1, kx = nlev[x]),
  posterior_sample.backdoor_params_cat(fit$sparse, 10**3, ess = 1, kx = nlev[x])
)

