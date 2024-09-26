

library(dplyr)
library(ggplot2)

dir_in <- "./inst/simulations/results/"
dir_out <- "./inst/simulations/plots/"
dir.create(dir_out)
files <- list.files("./inst/simulations/R", full.names = T)
sapply(files, source, echo = T)



res_from_file_to_df <- function(files) {
  imp <- lapply(files,
                function(f) do.call(c, readRDS(f)[c("par", "res")]))
  dfs <- lapply(imp, data.frame)
  df  <- do.call(rbind, dfs)
  names(df) <- gsub("^res.|^par.", "", names(df))
  names(df)[names(df) == "complexity"] <- "treedepth"

  for (v in c("n", "k", "treedepth")) {
    lev <- sort(unique(df[[v]]))
    lab <- paste0(v, "=", lev)
    df[[v]] <- factor(df[[v]], lev, lab)
  }

  df$N <- with(df, factor(N, sort(unique(N))))
  df$lstruct.epf <- with(df, interaction(local_struct, edgepf))
  return(df)
}



# n = 10
files <- list.files(dir_in, ".rds", full.names = T)
files <- files[!grepl("n8", files)]
files <- files[!grepl("N10000", files)]
df <- res_from_file_to_df(files)

# box-plots by N
x <- "N"
color <- "lstruct.epf"
facets <- "n+k~treedepth"
caption <- cat("The box-plots show the distribution over", paste0(range(df$r), collapse = "-"), "simulation runs.")
plots <- list()
for (y in c("edgep", "arp", "fpr", "tpr", "mse.unknown", "mse.known")) {
  plots[[y]] <- plot_boxplot(df, "N", y, facets, color = color, caption = caption)
}
lapply(names(plots), function(v) ggsave(paste0(dir_out, v, ".png"), plots[[v]], height = 5, width = 4))


title <- "Precision-recall of ancestor relation probabilities and edgeprobabilities. Closer to 1 is better. Averaged over all simulation runs."
file <- paste0(dir_out, "edgep_arp.tex")
keys  <- c("n", "k",  "treedepth")
values <- c("edgep", "arp")
method <- "lstruct.epf"
agg <- aggregate(df[values], df[c(keys, "N", method)], "mean")
values_from <- levels(agg[[method]])[unique(agg[[method]])]
agg %>%
  tidyr:::pivot_longer(values) %>%
  tidyr:::pivot_wider(names_from = method) %>%
  group_by(across(all_of(keys))) %>% #-> df
  arrange(N, .by_group = T) %>%
  df_to_tex(values_from = values_from,
            names_from = "name",
            caption = title,
            file = file)


title <- "MSE of intervention probabilities. With unknown and known structure. Closer to zero is better. Averaged over all simulation runs."
file <- paste0(dir_out, "mse.tex")
values <- c("mse.unknown", "mse.known")
keys   <- c("n", "k",  "treedepth")
method <- "lstruct.epf"
agg <- aggregate(df[values], df[c(keys, "N", method)], "mean")
values_from <- levels(agg[[method]])[unique(agg[[method]])]
agg %>%
  tidyr:::pivot_longer(values) %>%
  tidyr:::pivot_wider(names_from = method) %>%
  group_by(across(all_of(keys))) %>%
  df_to_tex(values_from = values_from,
            names_from = "name",
            caption = title,
            file = file)
