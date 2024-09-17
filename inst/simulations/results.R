

library(dplyr)
library(ggplot2)

dir_in <- "./inst/simulations/results/"
files <- list.files("./inst/simulations/R", full.names = T)
sapply(files, source, echo = T)



res_from_file_to_df <- function(files) {
  imp <- lapply(files,
                function(f) do.call(c, readRDS(f)[c("par", "res")]))
  dfs <- lapply(imp, data.frame)
  df  <- do.call(rbind, dfs)
  names(df) <- gsub("^res.|^par.", "", names(df))

  df$N <- with(df, factor(N, sort(unique(N))))
  df$complexity <- with(df, factor(complexity, sort(unique(complexity)), c("low", "high")))
  df$lstruct.epf <- with(df, interaction(local_struct, edgepf))
  return(df)
}



# n = 10
files <- list.files(dir_in, ".rds", full.names = T)
files <- files[!grepl("n8", files)]
df <- res_from_file_to_df(files)


x <- "N"
color <- "lstruct.epf"
facets <- "n+k~complexity"

plots <- list()
for (y in c("edgep", "arp", "fpr", "tpr", "mse.unknown", "mse.known")) {
  plots[[y]] <- plot_boxplot(df, "N", y, facets, color = color)
}
lapply(plots, print)
