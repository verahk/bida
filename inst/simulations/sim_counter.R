
branch  <- system("git branch --show-current", intern = TRUE)
paths   <- list(MCMC = paste0("./inst/simulations/", branch, "/MCMCchains/"),
                results = paste0("./inst/simulations/", branch, "/results/"))
files <- lapply(paths, list.files, pattern = ".rds")

tmp <- lapply(files, stringr::str_split, pattern = "_", simplify = TRUE)
dfs <- lapply(tmp, data.frame)
for (df in dfs) {
  print(as.data.frame(table(df[1:6])))
}

