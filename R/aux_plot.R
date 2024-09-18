
#' Title
#'
#' @param N (integer) sample size
#' @param r (integer) iteration number
#' @param parts (integer vector)
#' @param levels (list) levels of parent variables
#' @param k (integer) cardinality of outcome
#' @param verbose (bolean)
#'
#' @return matrix
#'
#' @examples
#'
#' levels <- list(0:1, 0:1)
#' parts  <- c(1, 1, 3, 4)
#' sim_run(N, r, parts, levels, k = 2, verbose = TRUE)
sim_run <- function(N, r, partname, n, k, verbose = F) {
  args <- as.list(environment())
  cat("\n", paste(names(args), args, sep = "="), "\n")

  tic  <- list(Sys.time())
  ess  <- 1
  nlev <- rep(k, n)

  # ground truth ----
  set.seed(r)
  partition <- sim_draw_partition(nlev, 1, partname)
  bn <- sim_draw_bn(nlev, 1, partition)
  parts <- bida:::get_parts(partition)
  tcpt <- matrix(bn[[1]]$prob, nrow = k**(n-1), ncol = k, byrow = T)
  pdo <- bida:::interv_probs_from_bn(bn, "exact")
  truth <- list(parts = parts,
                tcpt = tcpt,
                pdo = pdo)

  # simulate data ----
  set.seed(r+N)
  data <- bida:::sample_data_from_bn(bn, N)
  tic$data <- Sys.time()

  # compute counts and store in bdeu-object
  bdeu <- bida:::bida_bdeu(data, 1, 2:n, ess = ess, nlev = nlev)


  # apply optimization routines ----
  res <- list()
  fit_and_eval <- function(method, partition) {
    fit <- replace(bdeu, "partition", list(partition))
    res[[method]] <<- sim_eval_fit(fit, truth)
    tic[[method]] <<- Sys.time()
  }
  for (method in c("tree", "ptree")) {
    tmp <- bida:::optimize_bdeu(bdeu, method, regular = method == "ldag")
    fit_and_eval(method, tmp$partition)
  }

  if (!(k > 3 && n > 5)) {
    # there are none optimize_bdeu-routine with method = "pcart",
    # since this method is applied to the data not the frequency table
    tmp <- bida:::optimize_partition_from_data(data,
                                               1,
                                               2:n,
                                               ess = ess,
                                               nlev = nlev,
                                               method = "pcart",
                                               verbose = verbose)
    fit_and_eval("pcart", tmp$partition)
  }

  # add references - full & true ----
  fit_and_eval("full", as.list(seq_len(k**(n-1))-1))
  fit_and_eval("true", partition)

  # store res in matrix and add running times
  resmat <- cbind(do.call(rbind, res),
                  tic = c(as.double(diff(do.call(c, tic))[names(res)], units = "secs")))


  # attr(res, "tic") <-  diff(do.call(c, tic))
  # attr(res, "args") <- list(args)
  # return(res)
  # list output
  out <- c(args, list(method = names(res), resmat))
  attr(out, "tic") <- diff(do.call(c, tic))
  return(out)
}


sim_draw_partition <- function(nlev, j, partname) {
  # draw a partition of parent space of variable j
  nlev <- nlev[-j]
  q <- prod(nlev)
  n <- length(nlev)
  p <- switch(partname,
       "binary" = rep(seq.int(nlev[1]), each = q/nlev[1]),
       "vstruct"= {
         stride <- c(1, cumprod(nlev[-length(nlev)]))
         c(rep(1, sum(stride)), (sum(stride)+1):q)
         },
       "low" = bida:::rand_partition(nlev, 1, "tree", FALSE, nextsplitprob = function(x) 2/3),
       "high" = bida:::rand_partition(nlev, 1, "tree", FALSE, nextsplitprob = function(x) 1/3),
       "full" = seq_len(q))
  split(seq_along(p)-1, p)
}

sim_draw_bn <- function(nlev, j, partition) {
  # draw a bn with local structure for a small dag where
  n <- length(nlev)

  dag <- matrix(0, n, n)
  dag[-j, j] <- 1
  colnames(dag) <- rownames(dag) <- paste0("X", seq_len(n))
  partitions <- vector("list", n)
  partitions[[j]] <- partition
  dist <- bida:::rand_dist_cat(dag, 1, nlev, partitions)

  # replace marginal dist of parents with uniform
  for (i in seq_along(nlev)[-j]) dist[[i]][] <- 1/nlev[i]

  bn <- bida:::rand_bn(dag, "cat", alpha = 1, nlev, partitions)
  g <- bnlearn::empty.graph(colnames(dag))
  bnlearn::amat(g) <- dag
  bnlearn::custom.fit(g, dist)
}

sim_eval_fit <- function(fit, truth) {
  # compare fit to ground truth

  compute_paired_confusion_matrix <- function(x, y) {
    stopifnot(length(x) == length(y))
    n <- length(x)
    if (choose(n, 2) <= choose(2**8, 2)) {
      tmp <- combn(n, 2)
      pair_x <- x[tmp[1, ]] == x[tmp[2, ]]
      pair_y <- y[tmp[1, ]] == y[tmp[2, ]]

      out <- tabulate(pair_x + 2*pair_y +1, 4)
      names(out) <-  c("TN", "FP", "FN", "TP")
    } else {
      out <- setNames(rep(NA, 4), c("TN", "FP", "FN", "TP"))
      return(out)
      for (i in seq_len(n)) {
        cat("\niter:", i)
        indx <- seq_len(i)
        tmp <- 1 + (x[i] == x[-indx]) + 2*(y[i] == y[-indx])
        out <- out + tabulate(tmp, 4)
      }
    }
    out
  }

  parts <- bida:::get_parts(fit$partition)

  # compute point-estimates for the full CPT
  cpt <- bida:::posterior_mean(fit)[, parts]

  # compute point-estimates of the intervention probability, for each cause 2...n
  n   <- dim(truth$pdo)[1]
  pdo <- matrix(list(), n, n)
  for (i in seq(2, n)) {
    perm <- c(1, i, seq(2, n)[-(i-1)])
    tmp <- aperm(fit, perm)
    pdo[[i, 1]] <- bida:::backdoor_mean(tmp)
  }

  c(nparts = length(fit$partition),
    score  = bida:::score_bdeu(fit),
    sse_cpt    = sum((t(cpt)-truth$tcpt)**2),
    sse_ipt = sum((unlist(pdo[-1, 1])-unlist(truth$pdo[-1, 1]))**2),
    compute_paired_confusion_matrix(parts, truth$parts))
}


# test sim_run ----
if (FALSE) {
  sim_run(1000, 1, "binary", n = 3, k = 4)

  res <- lapply(1:2, function(r) sim_run(1000, r, "binary", n = 4+r, k = 4))
  res
}

# functions for summarizing and printing results ----
sim_res <- function(res, name_sim = NULL, dir_out = NULL) {
  # make plots, tables and save to file

  cat("Collect results in data frame\n")
  df <- sim_res_to_data_frame(res)

  # normalize
  vars_group = group_vars(df)
  vars_to_normalize <- c("nparts", "score", "sse_cpt", "sse_ipt")
  df_tmp <- df %>%
    group_by(across(all_of(c(vars_group[!vars_group == "method"], "r")))) %>%
    arrange(desc(method), .by_group = TRUE) %>%
    mutate(across(all_of(vars_to_normalize), ~ .x/first(.))) %>%
    # regrooup for correct aggregation below
    group_by(across(all_of(vars_group)))

  stopifnot(nrow(filter(df_tmp, method == "full" & nparts >1)) == 0)
  df <- df_tmp %>%
    filter(!method == "full")

  cat("Plot plots\n")
  plots <- list()
  plots[["nparts"]] <- df %>%
    sim_res_plot_boxplot(v = "nparts")
  plots[["score"]] <- df %>%
    sim_res_plot_boxplot(v = "score")
  plots[["sse_cpt"]] <- df %>%
    sim_res_plot_boxplot(v = "sse_cpt") +
    geom_abline(intercept = 1, slope = 0) +
    coord_cartesian(ylim = c(0, 3))
  plots[["sse_ipt"]] <- df %>%
    sim_res_plot_boxplot(v = "sse_ipt") +
    geom_abline(intercept = 1, slope = 0) +
    coord_cartesian(ylim = c(0, 3))
  plots[["dice"]] <- df %>%
    sim_res_plot_boxplot(v = "dice")
  plots[["pairs"]] <- df %>%
    sim_res_plot_pairs(average = TRUE)


  if (!is.null(dir_out)) {
    cat("Save plots\n")
    filenames <- paste0(name_sim, "_", names(plots), ".png")
    mapply(ggsave,
           filename = paste0(dir_out, filenames),
           plot = plots,
           height = 8, width = 6)


    # collect in markdown
    cat("Save results to .md file\n")
    file <- paste0(dir_out, name_sim, ".md")
    lines <- c("---",
               paste0("title: Simulation results:", name_sim),
               "output: html_document",
               "---",
               "\n",
               paste0(sprintf("\n![%s](%s)", filenames, filenames)))
    writeLines(lines, con = file)
    sink(file, TRUE)
    df %>%
      sim_res_table(c("score", "nparts", "sse_cpt", "sse_ipt", "tic")) %>%
      arrange(partname, name, n, k, N) %>%
      knitr:::kable(caption = "Averaged simulation results.", digits = 2)
    sink()

    rmarkdown::render(file, output_format = "html_document")
    browseURL(paste0(dir_out, name_sim, ".html"))

    # print latex tables
    cat("Print latex-tables\n")
    for (p in unique(df$partname)) {
      agg <- df  %>%
        filter(partname == p) %>%
        # compute mean and show as wide
        sim_res_table(c("score", "sse_cpt", "sse_ipt")) %>%
        group_by(n) %>%
        select(-all_of(c("partname", "k")))
      agg %>%  df_to_tex(names_from = "name",
                         values_from = c("pcart", "ptree", "tree", "true"),
                         file = paste0(dir_out, name_sim, "_", p, ".tex"),
                         caption = paste0(
                           "Simulation results. ",
                           sprintf("Cardinality: %s. Partition: %s.", unique(df$k), p)),
                         note = sprintf("All criteria relative to the full CPT and averaged over %s simulation runs.", max(df$r)),
                         addTimeStamp = T)
    }

    df %>%
      sim_res_table(c("tic")) %>%
      group_by(partname) %>%
      select(-true, -name, -k) %>%
      df_to_tex(names_from = "n",
                values_from = c("pcart", "ptree", "tree"),
                file = paste0(dir_out, name_sim, "_tic.tex"),
                group_keys_in_row = FALSE,
                caption = paste0(
                  "Runtimes in seconds. ",
                  sprintf("Cardinality: %s, Averaged over %s runs.", unique(df$k), max(df$r))),
                digits = 1)


  } else {
    return(plots)
  }
}


sim_res_to_data_frame <- function(res, groupvars = c("N", "n", "k", "partname", "method")) {
  # collect results in data frame
  vars_to_exclude <- c("verbose")
  df <- dplyr:::bind_rows(lapply(res, data.frame, row.names = NULL), .id = "id") %>%
    select(-any_of(vars_to_exclude))

  if ("TP" %in% names(df)) {
    df$dice = with(df, ifelse(TP == 0, 0, 2*TP/(2*TP+FP+FN)))
  }

  df %>%
    mutate(partname = factor(partname,
                             c("binary", "vstruct", "high", "low", "full")),
           method = factor(method, c("pcart", "ptree", "tree", "true", "full")),
           n = factor(n, n, paste0("n=", n)),
           k = factor(k, k, paste0("k=", k)),
           N = factor(N, sort(unique(N)))) %>%
    group_by(across(all_of(groupvars)))
}


sim_res_table <- function(df, vars) {
  df %>%
    select(vars) %>%
    summarize(across(vars, ~ mean(.x)), .groups = "keep") %>%
    tidyr::pivot_longer(vars) %>%
    tidyr::pivot_wider(names_from = "method") %>%
    group_by(partname, name) %>%
    select(partname, name, n, k, everything()) %>%
    arrange(n, k, N)
}


sim_res_prettify_plot <- function(xlab = NULL,
                                  ylab = NULL,
                                  caption = NULL) {
  list(theme_minimal(),
       theme(legend.position = "bottom",
             plot.caption = element_text(hjust = 0)),
       labs(caption = caption),
       ylab(ylab))
}

sim_res_plot_scatter <- function(df, x, y) {
  df %>%
    ggplot(aes(x = !! sym(x), y = !! sym(y),
               color = method, shape = as.factor(N))) +
    facet_grid(partname~n+k, scales = "free") +
    geom_point()
}

sim_res_plot_errorbar <- function(df, v) {
  df %>%
    tidyr::pivot_longer(all_of(v)) %>%
    group_by(name, .add = T) %>%
    summarize(high = quantile(value, .75),
              low = quantile(value, .25),
              value = median(value)) %>%
    ggplot(aes(N, y = value, ymin = low, ymax = high, color = method, linetype = name)) +
    facet_grid(partname~n+k, scales = "free") +
    geom_errorbar(width = .2) +
    geom_line(linewidth = .75) +
    sim_res_prettify_plot()
}
sim_res_plot_boxplot <- function(df, v) {
  caption <- sprintf("The box-plots shows the distribution over %s simulation runs.",
                     length(unique(df$r)))
  df %>%
    rename(value = !! sym(v)) %>%
    ggplot(aes(N, y = value, color = method, group = interaction(N, method))) +
    facet_grid(partname~n+k, scales = "free") +
    geom_boxplot() +
    #geom_jitter(width = 0.05, alpha = .7) +
    coord_cartesian(ylim = c(0, 1)) +
    sim_res_prettify_plot(ylab = v,
                          caption = caption)
}
sim_res_plot_pairs <- function(df, average = FALSE) {
  cols <-  c(TN = "grey",  FN = "blue", FP = "red",TP = "green")
  facets <- as.formula(n+k+partname+r ~ method)
  if (average) {
    df <- df %>%
      summarize(across(all_of(names(cols)), ~mean(.x)), .groups = "keep")
    facets <- as.formula(n+k+partname ~ method)
  }

  df %>%
    tidyr::pivot_longer(names(cols)) %>%
    mutate(name = factor(name, names(cols))) %>%
    #  summarize(sum(value)) %>% print(n = 1000)
    ggplot(aes(N, value, fill = name)) +
    facet_grid(facets, scales = "free") +
    geom_col() +
    scale_fill_manual(values = cols) +
    sim_res_prettify_plot()
}




