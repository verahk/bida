plot_box_plot_no_points <- function(grouped_df, x = "N", y = "value", fill = NULL, color = NULL, group = NULL, color_mean = NULL,  alpha = 1, width = .5) {

  xx <- sym(x)
  yy <- sym(y)
  if (!is.null(fill)) fill <- sym(fill)
  if (!is.null(color)) color <- sym(color)
  if (!is.null(color_mean) && !is.na(color_mean)) color_mean <- sym(color_mean)

  plot <- grouped_df %>%
    summarise_box_plot(yy) %>%
    ggplot(aes(!! xx, fill = !! fill, color = !! color, group = group))  + #fill = fill, color = color,
    geom_boxplot_no_points(alpha, width)

  if (!(!is.null(color_mean) && is.na(color_mean))) plot <- plot + geom_point(aes(y = mean, color = color_mean), stroke = .25, shape = 20)
  return(plot)
}

geom_boxplot_no_points <- function(width = .5, alpha = .75) {
  geom_boxplot(stat = "identity",
               # position = position_dodge2(padding = 0),
               width = width,
               alpha = alpha,
               aes(lower = p25,
                   upper = p75,
                   middle = p50,
                   ymin = pmax(p50-IQR*1.5, ymin),
                   ymax = pmin(p50+IQR*1.5, ymax)))
}

summarise_box_plot <- function(grouped_df, value = sym("value")) {
  grouped_df %>%
    summarize(mean = mean(value, na.rm = T),
              p25  = quantile(value, p = .25, na.rm = T),
              p50  = median(value, na.rm = T),
              p75  = quantile(value, p = .75, na.rm = T),
              IQR  = p75-p25,
              ymin = min(value),
              ymax = max(value),
              .groups = "keep")  %>%
    mutate(group = cur_group_id())
}


