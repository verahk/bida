
## theme
plot_adjust_theme <- function(axis.text.x = element_text(), xlab = "", ylab = "", col_values = NULL){
  x <- list(scale_fill_discrete(name = ""),
            scale_color_discrete(name = ""),
            theme_light(),
            theme(legend.position = "bottom",
                  axis.text.x = axis.text.x,
                  strip.text = element_text(colour = 'black'),
                  strip.background = element_rect(fill = 'white')),
            xlab(xlab),
            ylab(ylab))

  if (!is.null(col_values)) {
    x[1:2] <- list(scale_fill_manual(name = "", values = col_values),
                   scale_color_manual(name = "", values = col_values))
  }
  return(x)
}

