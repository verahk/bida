
theme_bida <- function(...){
  list(ggplot2::theme_light(),
       ggplot2::theme(legend.position = "bottom",
                      legend.title = ggplot2:::element_blank(),
                      strip.text = ggplot2:::element_text(colour = 'black'),
                      strip.background = ggplot2:::element_rect(fill = 'white'),
                      plot.caption = element_text(hjust = 0),
                      ...))
}

