

library(dplyr)
library(ggplot2)

dir_in <- "./inst/simulations/results/"
dir_out <- "./inst/simulations/plots/"
dir.create(dir_out)
files <- list.files("./inst/simulations/R", full.names = T)
sapply(files, source, echo = T)


# list files with results ----
files <- list.files(dir_in, ".rds", full.names = T)

res_from_file_to_df <- function(files, name) {
  imp <- lapply(files,
                function(f) do.call(c, readRDS(f)[c("par", name)]))
  dfs <- lapply(imp, data.frame)
  df  <- do.call(rbind, dfs)
  names(df) <- gsub(paste0(name, "."), "", names(df))
  names(df) <- gsub(paste0("par", "."), "", names(df))
  names(df)[names(df) == "complexity"] <- "csi"

  for (v in c("n", "k", "csi")) {
    lev <- sort(unique(df[[v]]))
    lab <- paste0(v, "=", lev)
    df[[v]] <- factor(df[[v]], lev, lab)
  }

  if (all(c("known", "unknown", "full") %in% names(df))) {
    df <- df %>%
      tidyr::pivot_longer(any_of(c("known", "unknown", "full", "arp"))) %>%
      mutate(slearn = ifelse(name == "known", "true parents", "unknown parents"),
             name   = ifelse(name == "full", "full CPT", "reduced CPT"),
             value  = sqrt(value))
  }

  df$N <- with(df, factor(N, sort(unique(N))))
  df$lstruct.epf <- with(df, interaction(local_struct, edgepf))
  return(df)
}



sim_results <- function(df, name, varsdoSave = FALSE) {
  keys <- c("init", "local_struct", "sample", "ess", "edgepf", "hardlimit",
            "N", "n", "k", "csi", "r", "lstruct.epf")
  vars <- names(df)[!names(df) %in% keys]
  # boxplots
  x <- "N"
  color <- "lstruct.epf"
  facets <- "n+k~csi"
  caption <- cat("The box-plots show the distribution over", paste0(range(df$r), collapse = "-"), "simulation runs.")
  plots <- list()
  for (y in vars) {
    plots[[y]] <- plot_boxplot(df, "N", y, facets, color = color, caption = caption)
  }
  if (doSave) {
    lapply(names(plots), function(v) ggsave(paste0(dir_out, name, "_", v, ".png"), plots[[v]], height = 5, width = 4))
  } else {
    plots
  }
}

stop()
res_from_file_to_df(files, "mse_pdo") %>%
  plot_boxplot("N", "value", "n+k~csi+slearn+name", color = "lstruct.epf", ylab = "RMSE, IPTs")
res_from_file_to_df(files, "rank") %>%
  plot_boxplot("N", "value", "n+k+csi~slearn+name", color = "lstruct.epf", ylab = "average precision")
res_from_file_to_df(files, "mse_tau") %>%
  plot_boxplot("N", "value", "n+k+csi~slearn+name", color = "lstruct.epf", ylab = "RMSE, causal effects")
res_from_file_to_df(files, "parents") %>%
  plot_boxplot("N", "value", "n+k+csi~slearn+name", color = "lstruct.epf", ylab = "size of conditioning set")
res_from_file_to_df(files, "parts") %>%
  plot_boxplot("N", "value", "n+k+csi~slearn+name", color = "lstruct.epf", ylab = "size of IPT (number of rows)")

df <- res_from_file_to_df(files, "mse_pdo")



# tabs ----
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
