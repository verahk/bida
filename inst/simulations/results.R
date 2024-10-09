

library(dplyr)
library(ggplot2)

do_print_tabs <- F

branch  <- system("git branch --show-current", intern = TRUE)
dir_in  <- paste0("./inst/simulations/", branch, "/results/")
dir_out <- paste0("./inst/simulations/", branch, "/output/")
dir.create(dir_out)
files <- list.files("./inst/simulations/R", full.names = T)
sapply(files, source, echo = T)

res_from_file_to_df <- function(files, name) {
  imp <- lapply(files,
                function(f) do.call(c, readRDS(f)[c("par", name)]))
  rename <- function(x) {
    namesx <- names(x)
    pos <- match(paste0(name, ".n"), namesx, 0L)
    if (pos>0) {
      # rename n at "positive edges"
      namesx[pos] <- paste0(namesx[pos], "pos")
      names(x) <- namesx
    }
    x
  }
  imp <- lapply(imp, rename)
  dfs <- lapply(imp, data.frame)
  df  <- do.call(rbind, dfs)
  names(df) <- gsub(paste0(name, "\\."), "", names(df))
  names(df) <- gsub(paste0("par", "\\."), "", names(df))



  if (all(c("known", "unknown", "full") %in% names(df))) {
    df <- df %>%
      tidyr::pivot_longer(any_of(c("known", "unknown", "full", "arp"))) %>%
      mutate(slearn = ifelse(name == "known", "true parents", "unknown parents"),
             name   = case_when(name == "full" ~ "full CPT",
                                name == "arp" ~ "ARP",
                                .default = "reduced CPT"))


  }

  for (v in c("n", "k", "N")) {
    df[[v]] <- with(df, factor(df[[v]], sort(unique(df[[v]]))))
  }


  df$maxdepth <- paste0("treedepth=", df$maxdepth)
  df$edgepf <- ifelse(df$edgepf == 2, "2", "logN")
  df$lstruct.epf <- with(df, interaction(local_struct, edgepf))

  # check that there are not more than 1 line per group
  group_vars <- c("n", "k",  "N", "maxdepth", "lstruct.epf", "slearn", "name")
  tmp <- group_by(df, across(any_of(c(group_vars, "r")))) %>% count(name = "count")
  stopifnot(all(tmp$count) == 1)

  # return grouped df - removing iter r
  df %>%
    select(group_vars, everything()) %>%
    group_by(across(any_of(group_vars)))
}


# Summarize ----

# Structure learning -----
if (do_print_tabs) {
files <- list.files(dir_in, ".rds", full.names = T)
files <- files[!grepl("none.*logN", files)]

values_from <- c("fpr", "tpr", "avgppv")
names_from  <- c("lstruct.epf", "maxdepth")
row_group_by <- c("n", "k")
pretty_names <- c(edgep = "edge", arp = "ancestor relation")

for (name in names(pretty_names)) {

  df <- res_from_file_to_df(files, name)%>%
        summarize(across(all_of(values_from), ~ mean(.x)),
                  nr = n(),
                  .groups = "keep")

  titles <- c(fpr = paste0("False positive rate of ", pretty_names[name], "s, averaged over all DAGs in the sample."),
              tpr = paste0("True positive rate of ", pretty_names[name], "s, averaged over all DAGs in the sample"),
              avgppv = paste0("Average precision using the posterior ", pretty_names[name], " probabilities."))
  titles <- paste(titles,
                   "Averaged over", round(mean(df$nr)), "simulation runs.")
  names(titles) <- c("fpr", "tpr", "avgppv")
  for (v in values_from) {
    file  <- paste0(dir_out, "slearn_", name, "_", v, ".tex")
    df %>%
      select(all_of(c(row_group_by, "N", v))) %>%
      group_by(across(all_of(row_group_by))) %>%
      arrange(N, .by_group = T) %>% #-> df
      df_to_tex(values_from = v,
                names_from = names_from,
                caption = titles[v],
                file = file,
                label = paste0("tab:",substr(file, nchar(dir_out)+1, nchar(file)-4)))

    file  <- paste0(dir_out, "slearn_", name, "_", v, "_n20.tex")
    df %>% filter(n == "20") %>%
      select(all_of(c(row_group_by, "N", v))) %>%
      group_by(across(all_of(row_group_by))) %>%
      arrange(N, .by_group = T) %>% #-> df
      df_to_tex(values_from = v,
                names_from = names_from,
                caption = titles[v],
                file = file,
                label = paste0("tab:",substr(file, nchar(dir_out)+1, nchar(file)-4)))
  }


}

# Parent set size ----
files <- list.files(dir_in, ".rds", full.names = T)
files <- files[grepl("n20", files)]
files <- files[!grepl("none.*logN", files)]

n_sim_runs <- format(mean(table(gsub("r[0-9]+", "", files))), digits = 1)

names_from  <- c("lstruct.epf", "maxdepth")
row_group_by <- c("n", "k")
title <- paste0("Average parent set size in the sampled DAGs.",
                " Averaged over ", n_sim_runs, " simulation runs.")
file  <- paste0(dir_out, "parent_set_size_n20.tex")


df <- res_from_file_to_df(files, "parents") %>%
  filter(!(name == "reduced CPT" | slearn == "known parents")) %>%
  summarize(value = mean(value), .groups = "drop") %>%
  select(-name, -slearn) %>%
  group_by(across(all_of(row_group_by))) %>%
  arrange(N, .by_group = T)

df %>%
  df_to_tex(values_from = "value",
            names_from = names_from,
            caption = title,
            file = file,
            label = paste0("tab:",substr(file, nchar(dir_out)+1, nchar(file)-4)))

# Run times ----

# Size of local backdoor distributions

}




path_report <- paste0(dir_out, "report.md")
write(Sys.time(), file = path_report)

# Causal-effect estimates ----
for (n in c(10, 20)) {
  for (k in c(2, 4)) {
  tag <- sprintf("n%s_k%s", n, k)
  files <- list.files(dir_in, ".rds", full.names = T)
  files <- files[grepl(tag, files)]
  if (length(files) < 30) next

  x <- "N"
  y <- "value"
  facets <- "n+k+maxdepth~slearn+name"
  fill <- "lstruct.epf"

  ylabs <- list("mse_pdo" = "intervention probs, MSE",
                "mse_tau" = "causal effects, MSE",
                "rank" = "positive causal effects, average precsision",
                "ranktop" = "top-causal effects, average precision",
                "parents" = "size of conditioning set, avg.",
                "parts" = "size of CPT, avg.")

  for (name in names(ylabs)) {
    plot <- res_from_file_to_df(files, name) %>%
      filter(!name == "full CPT") %>%
      plot_boxplot(x, y, facets, fill = fill, ylab = ylabs[[name]])
    if (length(dir_out) > 0 ) {
      filename <- paste0("box_plot_", sprintf("n%s_k%s", n, k), "_", name, ".png")
      ggsave(paste0(dir_out, filename), plot, height = 7, width = 6)

      new_line <- sprintf("![%s](%s)", filename, filename)
      write(new_line, file = path_report, append = TRUE)
    } else {
      print(plot)
    }
  }


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
}


rmarkdown::render(path_report, output_format = "html_document")
browseURL(gsub("md$", "html", path_report))

stop()
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
