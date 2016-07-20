plot_tax_group_contribution <- function(tax.gr.props, by.years = FALSE) {
  ## prepares and plots in ggplot2 the contributon of taxonomic groups to the 
  ## abundance/biomass (exressed as proportions) of the communities. 
  ## Arguments: tax.gr.props - data frame (in wide format) of taxonomic groups'
  ##              proportions; stations x taxonomic groups, with additional 
  ##              identifying factors (expected: stations and years; all others
  ##              are ignored by default).
  ##            by.years - should the final plot be faceted by year?
  ## Returns ggplot2 object.
  
  library(ggplot2)
  library(plyr)
  library(reshape2)
  library(viridis)
  
  if(by.years) {
    # summarize data by stations and years, if specified at input 
    summary.tax.props <- ddply(tax.gr.props, .(stations, years), colwise(mean, .cols = is.numeric))  
    
  } else {
    # summarize data by stations only
    summary.tax.props <- ddply(tax.gr.props, .(stations), colwise(mean, .cols = is.numeric))  
  }
  
  # reshape the data frame in long format
  tax.props.melted <- melt(summary.tax.props)  
  
  # plot proportions as stacked bars
  p <- ggplot(data = tax.props.melted, aes(x = stations, y = value, fill = variable)) + 
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(discrete = TRUE, name = "") +
    labs(x = "Station", y = "Proportion") + 
    theme_bw() + 
    # modify the size of x axis and legend labels 
    theme(axis.text.x = element_text(size = rel(1.3)),
          axis.text.y = element_text(size = rel(1.3)),
          axis.title.x = element_text(size = rel(1.2)),
          legend.text = element_text(size = rel(1.1)))
  
  # facet by year, if specified at input
  if(by.years){
    p <- p + facet_wrap(~years, nrow = 2) + 
             theme(strip.text.x = element_text(size = rel(1.3)))
  }
  
  return(p)
}
