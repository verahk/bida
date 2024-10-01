

library(dplyr)
library(ggplot2)

branch  <- system("git branch --show-current", intern = TRUE)
dir_in  <- paste0("./inst/simulations/", branch, "/results/")
dir_out <- paste0("./inst/simulations/", branch, "/")
dir.create(dir_out)
files <- list.files("./inst/simulations/R", full.names = T)
sapply(files, source, echo = T)


res_from_file_to_df <- function(files, name) {
  imp <- lapply(files,
                function(f) do.call(c, readRDS(f)[c("par", name)]))
  rename_local_struct <- function(x) {
    names <- names(x)
    pos <- match("par.regular", names, 0L)
    if (pos > 0 && x[pos] == 1) {
      x <- x[-pos]
      x["par.local_struct"] = paste0(x["par.local_struct"], "reg")
    }
    return(x)
  }
  imp <- lapply(imp, rename_local_struct)
  dfs <- lapply(imp, data.frame)
  df  <- do.call(rbind, dfs)
  names(df) <- gsub(paste0(name, "\\."), "", names(df))
  names(df) <- gsub(paste0("par", "\\."), "", names(df))
  names(df)[names(df) == "complexity"] <- "csi"

  df$csi <- factor(df$csi, c(0, .5, 1), c("high-complex", "low-complex", "no-CSI"))
  for (v in c("n", "k")) {
    lev <- sort(unique(df[[v]]))
    lab <- paste0(v, "=", lev)
    df[[v]] <- factor(df[[v]], lev, lab)
  }

  if (all(c("known", "unknown", "full") %in% names(df))) {
    df <- df %>%
      tidyr::pivot_longer(any_of(c("known", "unknown", "full", "arp"))) %>%
      mutate(slearn = ifelse(name == "known", "true parents", "unknown parents"),
             name   = case_when(name == "full" ~ "full CPT",
                                name == "arp" ~ "ARP",
                                .default = "reduced CPT"))


  }

  df$N <- with(df, factor(N, sort(unique(N))))
  df$lstruct.epf <- with(df, interaction(local_struct, edgepf))

  return(df)
}

path_report <- paste0(dir_out, "report.md")
write(Sys.time(), file = path_report)
for (k in c(2, 4, 8)) {
  files <- list.files(dir_in, ".rds", full.names = T)
  files <- files[grepl(paste0("k", k), files)]

  x <- "N"
  y <- "value"
  facets <- "n+k+csi~slearn+name"
  color  <- "lstruct.epf"

  ylabs <- list("mse_pdo" = "intervention probs, MSE",
                "mse_tau" = "causal effects, MSE",
                "rank" = "positive causal effects, average precsision",
                "ranktop" = "top-causal effects, average precision",
                "parents" = "size of conditioning set, avg.",
                "parts" = "size of CPT, avg.")

  for (name in names(ylabs)) {
    plot <- res_from_file_to_df(files, name) %>%
      plot_boxplot(x, y, facets, color = color, ylab = ylabs[[name]])
    if (length(dir_out) > 0 ) {
      filename <- paste0("box_plot_k", k, "_", name, ".png")
      ggsave(paste0(dir_out, filename), plot, height = 7, width = 6)

      new_line <- sprintf("![%s](%s)", filename, filename)
      write(new_line, file = path_report, append = TRUE)
    } else {
      print(plot)
    }
  }

  rmarkdown::render(path_report, output_format = "html_document")
  browseURL(gsub("md$", "html", path_report))

  if (FALSE) {
    gr_vars <- c("init", "local_struct", "sample", "ess", "edgepf", "hardlimit",
                 "n", "k", "csi", "lstruct.epf", "slearn", "name", "N")

    df <- res_from_file_to_df(files, "size") %>%
      select(-r, -N, -lstruct.epf) %>%
      group_by(across(any_of(gr_vars))) %>%
      slice(1)

    # tab: average set of conditioning set
    title <- "Local distribution P(Y|X, Z). Size of conditioning set and number of rows. Averaged over all DAGs and nodes."
    file <- paste0(dir_out, "tab_cpt_size_k", k, "_", name, ".tex")
    dfs <- list(vars = res_from_file_to_df(files, "parents"),
                parts = res_from_file_to_df(files, "parts"))
    # aggregate
    agg <- dplyr::bind_rows(dfs, .id = "variable") %>%
      select(-ess, -init, -sample, -local_struct, -edgepf) %>%
      group_by(across(any_of(c(gr_vars, "variable")))) %>%
      summarize(value = mean(value), .groups = "keep") %>%
      arrange()

    # print to tex
    df <- agg %>%
      tidyr::pivot_wider(names_from = "variable") %>%
      group_by(n, k, csi, slearn, name)

    df %>%
      df_to_tex(values_from = unique(agg$variable),
                names_from = "lstruct.epf",
                caption = title,
                file = file)
  }
}




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
