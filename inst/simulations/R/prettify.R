
prettify_plot <- function(xlab = NULL,
                          ylab = NULL,
                          caption = NULL) {
  list(theme_minimal(),
       theme(legend.position = "bottom",
             plot.caption = element_text(hjust = 0)),
       labs(caption = caption),
       ylab(ylab))
}

