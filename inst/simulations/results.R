


dir_in <- "./inst/simulations/results/"
files <- list.files(dir_in, ".rds", full.names = T)
sim_res_to_df <- function(files) {
  imp <- lapply(files,
                function(f) do.call(c, readRDS(f)[c("par", "res")]))
  dfs <- lapply(imp, data.frame)

  df  <- do.call(rbind, dfs)
  names(df) <- gsub("^res.|^par.", "", names(df))
}

