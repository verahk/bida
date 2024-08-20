


# Test that `part` and `ldag` optimal partition routines return the same partition
levels <- rep(list(0:1), 4)
nlev <- lengths(levels)
ess <- 1

for (i in 1:100) {
  #cat("\niter:", i, "\n")
  counts <- matrix(rgamma(2**(length(levels)+1), shape = 1), ncol = 2)
  fits <- list()

  fits$part <- optimize_partition(counts, levels, ess, "part", regular = T, F)
  fits$ldag <- optimize_partition(counts, levels, ess, "ldag", regular = T, F)

  partitions <- sapply(lapply(fits, "[[", "partition"), get_parts)


  # stop and print output if part and ldag do not coincide
  if (!all(partitions[, "part"] == partitions[, "ldag"])) {
    cbind(expand.grid(levels), partitions)

    cat("\npart:")
    fits$part <- optimize_partition(counts, levels, ess, "part", regular = T, T)
    cat("\nldag:")
    fits$ldag <- optimize_partition(counts, levels, ess, "ldag", regular = T, T)

    is_CSI_consistent(fits$part$partition, nlev)
    is_regular(fits$part$partition, nlev)

    is_CSI_consistent(fits$ldag$partition, nlev)
    is_regular(fits$ldag$partition, nlev)
    break
  }
}
