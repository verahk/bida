

# args ----
args <- as.list(commandArgs(trailingOnly = TRUE))
if (length(args) == 0) {
  args <- list("bida", 0, 6)
} else if (length(args) < 3) stop()
names(args) <- c("what", "test_row", "nClusters")
args[-1] <- lapply(args[-1], as.numeric)
print(args)
vapply(args, class, character(1))

# paths ----
branch  <- system("git branch --show-current", intern = TRUE)
subdir  <- switch(args$what, "MCMC" = "MCMCchains", "bida" = "results")
outdir <- paste0("./inst/simulations/", branch, "/", subdir, "/")
subdir  <- switch(args$what, "bida" = "MCMCchains/")
indir <- paste0("./inst/simulations/", branch, "/", subdir, "/")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

branch  <- system("git branch --show-current", intern = TRUE)
paths   <- list(MCMC = paste0("./inst/simulations/", branch, "/MCMCchains/"),
                results = paste0("./inst/simulations/", branch, "/results/"))
files <- lapply(paths, list.files, pattern = ".rds")

tmp <- lapply(files, stringr::str_split, pattern = "_", simplify = TRUE)
dfs <- lapply(tmp, data.frame)
for (df in dfs) {
  print(as.data.frame(table(df[1:8])))
}

