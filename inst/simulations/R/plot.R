
if (FALSE) {
  library(extrafont)
  loadfonts()
  font_import(pattern = "lmroman*")
  # link www.fontsquirrel.com/fonts/latin-modern-roman

  # execute once to add fonts:
  # font_import(pattern = "lmroman*")

  # example
  df <- data.frame(x = 1:10, y = 1:10)
  ggplot(df, aes(x, y)) +
    theme(text = element_text(size=10, family="LM Roman 10"))

  plot_init(df, x, y)
}

plot_init <- function(df, x, y,
                      color = "",
                      fill = "",
                      shape = "",
                      linetype = "",
                      ylab = y,
                      xlab = x,
                      title = "",
                      subtitle = "",
                      caption = "") {

  ggplot(df, aes(x = !! sym(x),
                 y =  !! sym(y),
                 color = !! sym(color),
                 fill = !! sym(fill),
                 linetype = !! sym(linetype),
                 shape = !! sym(shape))) +
    labs(title = title,
         subtitle = subtitle,
         caption = caption) +
    ylab(ylab) +
    xlab(xlab) +
    plot_prettify()
}

plot_prettify <- function() {
  list(theme_minimal(),
       theme(legend.position = "bottom",
             plot.caption = element_text(hjust = 0)))
}

plot_add_facet_grid <- function(string, ...) {
  if (! (is.null(string) || string == "")) {
    list(facet_grid(as.formula(string), ...))
  }
}


plot_scatter <- function(df, x, y, facets = "local_struct~.", ...) {
  plot_init(df, x, y, facets = facets) +
    geom_point()
}

plot_errorbar <- function(df, x, y, facets = "", ...) {
  df %>%
    tidyr::pivot_longer(all_of(y)) %>%
    group_by(name, .add = T) %>%
    summarize(high = quantile(value, .75),
              low = quantile(value, .25),
              value = median(value)) %>%
    plot_init(x, y = "value", ...) +
      plot_add_facet_grid(facets, scale = "free") +
      geom_errorbar(aes(ymin = low, ymax = high), width = .2) +
      geom_line(linewidth = .75)
}

plot_boxplot <- function(df, x, y, facets = "", ...) {
    plot_init(df, x, y, ...) +
    plot_add_facet_grid(facets, scales = "free") +
    geom_boxplot()
    #geom_jitter(width = 0.05, alpha = .7)
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




1
