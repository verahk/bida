

sim_load_params <- function(which, what, iterStart, iterStop) {
  if (what == "MCMC") {

    par <- list(local_struct = c("ptree", "none"),
                init = c("pcskel"),
                sample = "order",
                ess = 1,
                edgepf = c("2", "logN"),
                hardlimit = 4,
                N = c(300, 1000, 3000),
                r = seq.int(iterStart, iterStop))

    if (which == "syntethic") {

      # add params controlling random-cpt generation
      par$n <- c(10, 20)
      par$k <- c(2, 4)
      par$maxdepth <- c(0, .5, 1)

      params_to_filename <<- function(par) {
        tmp <- sprintf("n%s_k%s_depth%s_%s_%s_%s_ess%s_epf%s_N%s_r%02.0f.rds",
                       par$n, par$k, par$maxdepth*100, par$init, par$local_struct, par$sample, par$ess, par$edgepf, par$N, par$r)
        stopifnot(length(tmp) == 1)  # fails if any argument is NULL
        tmp
      }
    } else {
      stop()
      par$bnname <- c("asia")

      params_to_filename <<- function(par) {
        tmp <- sprintf("%s_%s_%s_%s_ess%s_epf%s_N%s_r%02.0f.rds",
                       par$bnname, par$init, par$local_struct, par$sample, par$ess, par$edgepf, par$N, par$r)
        stopifnot(length(tmp) == 1)  # fails if any argument is NULL
        tmp
      }
    }
    expand.grid(par, stringsAsFactors = FALSE)
  }  else if (what == "bida") {
    stop()
    files  <- list.files(indir, ".rds")
    cbind(files)
  }
}

