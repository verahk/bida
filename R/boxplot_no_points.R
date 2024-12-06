#' Boxplot with no outliers
#'
#' @param df (data.frame)
#' @param x,y,fill,color,group (character) name of variables to corresponding arguments of [ggplot2::aes()]
#' @param color_mean (character)
#' @param alpha (numeric)
#' @param width (numeric)
#'
#' @return a ggplot2-object
#' @keywords internal
boxplot_no_points <- function(df,
                              x = "N",
                              y = "value",
                              fill = NULL,
                              color = NULL,
                              group = NULL,
                              color_mean = NULL,
                              alpha = 1,
                              width = .5) {

  if (!is.null(fill)) fill <- sym(fill)
  if (!is.null(color)) color <- sym(color)
  if (!is.null(color_mean) && !is.na(color_mean)) color_mean <- sym(color_mean)

  grouped_df <- dplyr::group_by(df, across(all_of(c(x, fill, color, group))))
  plot <- summarise_box_plot(grouped_df, y) %>% # summary stats
    ggplot(aes(!! sym(x), fill = !! fill, color = !! color, group = group)) +
    geom_boxplot(stat = "identity",
                 # position = position_dodge2(padding = 0),
                 width = width,
                 alpha = alpha,
                 aes(lower = p25,
                     upper = p75,
                     middle = p50,
                     ymin = pmax(p50-(p75-p25)*1.5, min),
                     ymax = pmin(p50+(p75-p25)*1.5, max)))

  if (!is.null(color_mean) && !is.na(color_mean)) plot <- plot + geom_point(aes(y = mean, color = color_mean), stroke = .25, shape = 20)
  return(plot)
}

summarystat_by_group <- function(grouped_df, varname = "value") {
  dplyr::summarize(grouped_df,
                    mean = mean(!! sym(varname), na.rm = T),
                    p25  = quantile(!! sym(varname), p = .25, na.rm = T),
                    p50  = median(!! sym(varname), na.rm = T),
                    p75  = quantile(!! sym(varname), p = .75, na.rm = T),
                    min = min(!! sym(varname)),
                    max = max(!! sym(varname)),
                    group = dplyr::cur_group_id(),
                    .groups = "keep")
}

