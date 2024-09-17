

plot_scatter <- function(df, x, y) {
  df %>%
    ggplot(aes(x = !! sym(x), y = !! sym(y),
               color = method, shape = as.factor(N))) +
    facet_grid(partname~n+k, scales = "free") +
    geom_point()
}

plot_errorbar <- function(df, v) {
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
plot_boxplot <- function(df, v) {
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
plot_pairs <- function(df, average = FALSE) {
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




