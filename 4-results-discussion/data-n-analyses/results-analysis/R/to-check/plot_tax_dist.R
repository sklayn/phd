plot_tax_distinctness <- function(tax.dist.data, x, y, group) {
  # given a data frame, plots the selected taxonomic diversity indices
  
  library(ggplot2)
  
  ggplot(tax.dist.data, aes(x = get(x), y = get(y))) +
    geom_point(colour = get(group)) +                   # color points by group
    labs(x = x, y = y, colour = group) +                # proper axis and legend labels
    theme_bw()                      

}

